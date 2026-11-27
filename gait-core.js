/*!
 * pelviswalk-core — 스마트폰 골반(천골) 착용 IMU 기반 10m 보행 분석 엔진
 * License: MIT  (c) 건운사팀 (최예겸, 이상석, 김채연, 방아인)
 *
 * 입력: 원시 센서 샘플 { t(ms), ax, ay, az (가속도, 중력 포함, m/s²), gx, gy, gz (각속도, deg/s, 기기 x·y·z축) }
 * 출력: 시공간 지표, 좌우 지표, 골반 3축 ROM, 추정 수직 지면반발력, 대칭성·규칙성, 품질 검사
 * 본 결과는 참고용이며 의학적 진단을 대체하지 않습니다.
 */
(function (root) {
  'use strict';
  const FS = 100;           // 리샘플링 주파수 (Hz)
  const DT = 1 / FS;

  // ---------- 기본 수학 ----------
  const mean = a => a.reduce((s, v) => s + v, 0) / (a.length || 1);
  const sd = a => { const m = mean(a); return Math.sqrt(mean(a.map(v => (v - m) ** 2))); };
  const median = a => { const s = [...a].sort((x, y) => x - y); const n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : NaN; };
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm = a => Math.hypot(a[0], a[1], a[2]);
  const scale = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const unit = a => scale(a, 1 / (norm(a) || 1));
  const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

  // 선형 보간 리샘플링
  function resample(t, ys, fs) {
    const t0 = t[0], t1 = t[t.length - 1];
    const n = Math.floor((t1 - t0) * fs / 1000) + 1;
    const tt = new Array(n), out = ys.map(() => new Array(n));
    let j = 0;
    for (let i = 0; i < n; i++) {
      const ti = t0 + i * 1000 / fs; tt[i] = ti;
      while (j < t.length - 2 && t[j + 1] < ti) j++;
      const span = t[j + 1] - t[j] || 1, w = clamp((ti - t[j]) / span, 0, 1);
      for (let k = 0; k < ys.length; k++) out[k][i] = ys[k][j] * (1 - w) + ys[k][j + 1] * w;
    }
    return { tt, out };
  }

  // 2차 Butterworth 저역통과 (양방향 = 4차, 위상 지연 없음)
  function butterLP(x, fc, fs) {
    const k = Math.tan(Math.PI * fc / fs), q = Math.SQRT2;
    const n = 1 / (1 + k * q + k * k);
    const b0 = k * k * n, b1 = 2 * b0, b2 = b0, a1 = 2 * (k * k - 1) * n, a2 = (1 - k * q + k * k) * n;
    const pass = s => {
      const y = new Array(s.length); let x1 = s[0], x2 = s[0], y1 = s[0], y2 = s[0];
      for (let i = 0; i < s.length; i++) {
        const v = b0 * s[i] + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1; x1 = s[i]; y2 = y1; y1 = v; y[i] = v;
      }
      return y;
    };
    return pass(pass(x).reverse()).reverse();
  }

  // 가우시안 평활 (sigma: 초)
  function gauss(x, sigma) {
    const s = sigma * FS, r = Math.ceil(3 * s), w = [];
    for (let i = -r; i <= r; i++) w.push(Math.exp(-0.5 * (i / s) ** 2));
    const ws = w.reduce((a, b) => a + b, 0);
    return x.map((_, i) => { let acc = 0; for (let j = -r; j <= r; j++) acc += x[clamp(i + j, 0, x.length - 1)] * w[j + r]; return acc / ws; });
  }
  const deriv = x => x.map((_, i) => (x[Math.min(i + 1, x.length - 1)] - x[Math.max(i - 1, 0)]) / (2 * DT));

  function movingSD(x, win) {
    const h = Math.floor(win / 2);
    return x.map((_, i) => sd(x.slice(Math.max(0, i - h), Math.min(x.length, i + h + 1))));
  }

  function findPeaks(x, minDist, minH, lo = 0, hi = x.length - 1) {
    const cand = [];
    for (let i = Math.max(1, lo); i < Math.min(x.length - 1, hi); i++)
      if (x[i] > x[i - 1] && x[i] >= x[i + 1] && x[i] > minH) cand.push(i);
    cand.sort((a, b) => x[b] - x[a]);
    const keep = [];
    for (const c of cand) if (keep.every(k => Math.abs(k - c) >= minDist)) keep.push(c);
    return keep.sort((a, b) => a - b);
  }

  function autocorr(x, lag) {
    const m = mean(x), n = x.length - lag; let num = 0, den = 0;
    for (let i = 0; i < x.length; i++) den += (x[i] - m) ** 2;
    for (let i = 0; i < n; i++) num += (x[i] - m) * (x[i + lag] - m);
    return den ? num / den : 0;
  }

  // 구간 내 선형 추세 제거 후 적분 → 각도/변위
  function integrateDetrend(x, a, b) {
    const seg = x.slice(a, b + 1), m = mean(seg);
    let acc = 0; const out = seg.map(v => (acc += (v - m) * DT));
    const n = out.length - 1, slope = n ? (out[n] - out[0]) / n : 0;
    return out.map((v, i) => v - out[0] - slope * i);
  }
  const normCurve = (seg, n = 101) => Array.from({ length: n }, (_, i) => {
    const p = i * (seg.length - 1) / (n - 1), j = Math.floor(p), w = p - j;
    return seg[j] * (1 - w) + (seg[Math.min(j + 1, seg.length - 1)]) * w;
  });
  const avgCurves = cs => cs.length ? cs[0].map((_, i) => mean(cs.map(c => c[i]))) : [];
  const si = (l, r) => (l + r) ? Math.abs(l - r) / (0.5 * (l + r)) * 100 : 0;

  // ---------- 메인 분석 ----------
  function analyze(raw, subject = {}, opts = {}) {
    const S = raw.samples;
    if (!S || S.length < 50) return { ok: false, error: '데이터가 너무 짧습니다.' };
    const distance = opts.distance || 10;
    const t = S.map(s => s.t);
    const dts = t.slice(1).map((v, i) => v - t[i]);
    const mdt = median(dts);
    const fsRaw = 1000 / mdt;
    const missing = dts.filter(d => d > 2.5 * mdt).reduce((s, d) => s + d / mdt - 1, 0) / S.length * 100;
    const hasGyro = S.some(s => s.gx || s.gy || s.gz);

    const { out } = resample(t, ['ax', 'ay', 'az', 'gx', 'gy', 'gz'].map(k => S.map(s => s[k] || 0)), FS);
    const [ax, ay, az, gx, gy, gz] = out.map(c => butterLP(c, 20, FS));
    const N = ax.length;

    // 중력 방향: 정지 보정 값 → 없으면 처음 0.8초 평균
    let g = raw.calib && raw.calib.g;
    if (!g) { const k = Math.min(N, 80); g = [mean(ax.slice(0, k)), mean(ay.slice(0, k)), mean(az.slice(0, k))]; }
    const G = norm(g), gU = unit(g);
    // 진행 방향: 기기 z축(화면 수직)을 수평면에 투영. 화면이 바깥(뒤)을 향하면 -z가 전방
    let zh = sub([0, 0, 1], scale(gU, gU[2]));
    if (norm(zh) < 0.3) zh = sub([0, 1, 0], scale(gU, gU[1]));
    const apU = scale(unit(zh), -1), mlU = cross(gU, apU);

    const aV = [], bw = [], wV = [], wAP = [], wML = [];
    for (let i = 0; i < N; i++) {
      const a = [ax[i], ay[i], az[i]], w = [gx[i], gy[i], gz[i]];
      const v = dot(a, gU); aV.push(v - G); bw.push(v / G);
      wV.push(dot(w, gU)); wAP.push(dot(w, apU)); wML.push(dot(w, mlU));
    }

    // 보행 구간 검출
    const msd = movingSD(aV, 50);
    const th = Math.max(0.45, 0.3 * Math.max(...msd));
    let seg0 = msd.findIndex(v => v > th), seg1 = N - 1 - [...msd].reverse().findIndex(v => v > th);
    if (seg0 < 0) return { ok: false, error: '보행이 감지되지 않았습니다.', quality: { fsRaw, missing } };

    // 걸음 주기 추정 (자기상관)
    const segA = aV.slice(seg0, seg1 + 1);
    let bestLag = 55, bestR = -1;
    for (let lag = 30; lag <= 90 && lag < segA.length / 2; lag++) { const r = autocorr(segA, lag); if (r > bestR) { bestR = r; bestLag = lag; } }
    const stepEst = bestLag * DT;

    // 걸음(IC) 검출: 평활 수직가속도 피크 → 직전 최대 기울기 지점을 초기접지(IC)로
    const sm = gauss(aV, 0.04), dsm = deriv(gauss(aV, 0.02));
    const segSm = sm.slice(seg0, seg1 + 1);
    const pk = findPeaks(sm, Math.round(bestLag * 0.6), mean(segSm) + 0.35 * sd(segSm), seg0, seg1);
    const steps = pk.map(p => {
      // 최대 기울기 지점에서 거슬러 올라가 상승 시작점(onset)을 초기접지로
      let ms = p, best = -Infinity;
      for (let i = Math.max(0, p - 20); i <= p; i++) if (dsm[i] > best) { best = dsm[i]; ms = i; }
      let ic = ms;
      while (ic > Math.max(0, ms - 8) && dsm[ic - 1] > 0.4 * best) ic--;
      return { ic, pk: p };
    });
    const nSteps = steps.length;
    if (nSteps < 4) return { ok: false, error: '걸음이 충분히 검출되지 않았습니다.', quality: { fsRaw, missing, steps: nSteps } };

    // 좌우 구분: IC 직후 골반 수직축 회전 방향. 첫 걸음 = 오른발(출발 안내) 기준으로 부호 매핑
    const signs = steps.map(s => Math.sign(mean(wV.slice(s.ic, Math.min(N, s.ic + 12)))) || 1);
    let alt = 0; for (let k = 1; k < nSteps; k++) if (signs[k] !== signs[k - 1]) alt++;
    const altRatio = alt / (nSteps - 1);
    // 기본: '오른발로 출발' 안내에 따라 첫 걸음 = 오른발, 이후 교대.
    // 자이로 회전 방향이 교대 패턴과 일관된지(altRatio)는 품질 정보로만 제공
    let side = steps.map((_, k) => (k % 2 === 0 ? 'R' : 'L'));
    if (opts.flipLR) side = side.map(s => (s === 'R' ? 'L' : 'R'));
    steps.forEach((s, k) => { s.side = side[k]; s.t = s.ic * DT; });

    // 안정 구간: 앞뒤 2보 제외 (걸음이 적으면 1보)
    const trim = nSteps >= 9 ? 2 : nSteps >= 6 ? 1 : 0;
    const inSteady = k => k >= trim && k < nSteps - trim;

    // 시간 지표
    const stepT = [];
    for (let k = 0; k < nSteps - 1; k++) stepT.push({ k, v: (steps[k + 1].ic - steps[k].ic) * DT, side: steps[k + 1].side });
    const steadyStepT = stepT.filter(s => inSteady(s.k + 1) && inSteady(s.k)).map(s => s.v);
    const useStepT = steadyStepT.length >= 2 ? steadyStepT : stepT.map(s => s.v);
    const cadence = 60 / mean(useStepT);
    const duration = (steps[nSteps - 1].ic - steps[0].ic) * DT;
    const speed = distance / duration;

    // 접촉시간: 같은 발 IC → (반대발 IC 직후의) 이지(FC) ≈ 반대발 하중 피크
    const contact = [], dsArr = [];
    for (let k = 0; k < nSteps - 1; k++) {
      const fc = steps[k + 1].pk;
      const tc = (fc - steps[k].ic) * DT;
      const stride = k + 2 < nSteps ? (steps[k + 2].ic - steps[k].ic) * DT : null;
      if (inSteady(k) && inSteady(k + 1)) {
        contact.push({ side: steps[k].side, v: tc });
        if (stride) dsArr.push(((steps[k].pk - steps[k].ic) + (steps[k + 1].pk - steps[k + 1].ic)) * DT / stride * 100);
      }
    }
    const cL = mean(contact.filter(c => c.side === 'L').map(c => c.v)), cR = mean(contact.filter(c => c.side === 'R').map(c => c.v));
    const strideTs = [];
    for (let k = 0; k + 2 < nSteps; k++) if (inSteady(k) && inSteady(k + 2)) strideTs.push((steps[k + 2].ic - steps[k].ic) * DT);
    const strideT = strideTs.length ? mean(strideTs) : 2 * mean(useStepT);
    const gctBalance = cL / (cL + cR) * 100;
    const stancePct = mean([cL, cR]) / strideT * 100;

    // 보폭: 역진자 모델 + 총합 10m 보정
    const legL = subject.legLength || (subject.height ? 0.53 * subject.height / 100 : 0.9);
    const rawSL = [];
    for (let k = 0; k < nSteps - 1; k++) {
      const a = steps[k].ic, b = steps[k + 1].ic;
      const vel = integrateDetrend(aV, a, b);
      let p = 0; const pos = vel.map(v => (p += v * DT));
      const n = pos.length - 1, sl = n ? (pos[n] - pos[0]) / n : 0;
      const pd = pos.map((v, i) => v - pos[0] - sl * i);
      const h = clamp(Math.max(...pd) - Math.min(...pd), 0.002, 0.15);
      rawSL.push({ side: steps[k + 1].side, v: 2 * Math.sqrt(Math.max(0, 2 * legL * h - h * h)), k });
    }
    const K = distance / rawSL.reduce((s, x) => s + x.v, 0);
    rawSL.forEach(x => (x.v *= K));
    const slSteady = rawSL.filter(x => inSteady(x.k) && inSteady(x.k + 1));
    const slUse = slSteady.length >= 2 ? slSteady : rawSL;
    const slL = mean(slUse.filter(x => x.side === 'L').map(x => x.v)), slR = mean(slUse.filter(x => x.side === 'R').map(x => x.v));

    // 골반 각도 (자이로 적분, 활보별 추세 제거)
    const pelvis = { available: hasGyro };
    const strideIdx = [];
    for (let k = 0; k + 2 < nSteps; k++) if (inSteady(k) && inSteady(k + 2) && steps[k].side === 'R') strideIdx.push([steps[k].ic, steps[k + 2].ic]);
    if (!strideIdx.length) for (let k = 0; k + 2 < nSteps; k += 2) strideIdx.push([steps[k].ic, steps[k + 2].ic]);
    if (hasGyro) {
      for (const [key, sig] of [['rotation', wV], ['obliquity', wAP], ['tilt', wML]]) {
        const curves = [], roms = [];
        for (const [a, b] of strideIdx) {
          const ang = integrateDetrend(sig, a, b);
          roms.push(Math.max(...ang) - Math.min(...ang));
          curves.push(normCurve(ang));
        }
        pelvis[key] = { rom: mean(roms), curve: avgCurves(curves) };
      }
    }

    // 추정 수직 지면반발력 (BW)
    const bwS = gauss(bw, 0.015);
    const peaks = steps.filter((_, k) => inSteady(k)).map(s => {
      const w = Math.round(0.35 * mean(useStepT) * FS);
      let mi = s.ic, mv = -Infinity;
      for (let i = s.ic; i < Math.min(N, s.ic + w); i++) if (bwS[i] > mv) { mv = bwS[i]; mi = i; }
      return { side: s.side, peak: mv, lr: mv / Math.max(0.02, (mi - s.ic) * DT) };
    });
    const grfCurve = avgCurves(strideIdx.map(([a, b]) => normCurve(bwS.slice(a, b + 1))));
    const pkL = mean(peaks.filter(p => p.side === 'L').map(p => p.peak)), pkR = mean(peaks.filter(p => p.side === 'R').map(p => p.peak));
    const lrMean = mean(peaks.map(p => p.lr));

    // 대칭성·규칙성·변동성
    const stL = mean(stepT.filter(s => s.side === 'L' && inSteady(s.k + 1)).map(s => s.v)), stR = mean(stepT.filter(s => s.side === 'R' && inSteady(s.k + 1)).map(s => s.v));
    const Ad1 = autocorr(segA, bestLag), Ad2 = autocorr(segA, Math.min(segA.length - 2, 2 * bestLag));
    const cv = strideTs.length >= 2 ? sd(strideTs) / mean(strideTs) * 100 : sd(useStepT) / mean(useStepT) * 100;

    const quality = {
      fsRaw, missing, steps: nSteps, lrConsistency: hasGyro ? altRatio : null,
      ok: nSteps >= 6 && Math.round(fsRaw) >= 50 && missing <= 5,
      issues: [
        nSteps < 6 && `걸음 수 부족 (${nSteps}보 < 6보)`,
        Math.round(fsRaw) < 50 && `샘플링 속도 낮음 (${fsRaw.toFixed(0)} Hz < 50 Hz)`,
        missing > 5 && `결측률 높음 (${missing.toFixed(1)}% > 5%)`,
        !hasGyro && '자이로스코프 미지원: 골반 각도 생략',
      ].filter(Boolean),
    };

    return {
      ok: true, quality,
      general: { duration, steps: nSteps, cadence, speed, stepTime: mean(useStepT) },
      gct: { left: cL, right: cR, balance: gctBalance },
      stride: { left: slL, right: slR, strideLength: slL + slR, strideTime: strideT, stancePct, swingPct: 100 - stancePct, doubleSupportPct: mean(dsArr) },
      pelvis,
      grf: { left: pkL, right: pkR, loadingRate: lrMean, curve: grfCurve },
      symmetry: { stepTime: si(stL, stR), contact: si(cL, cR), stepLength: si(slL, slR), grf: si(pkL, pkR), Ad1, Ad2, regularityRatio: Ad1 / (Ad2 || 1), cv },
      trace: { t: aV.map((_, i) => i * DT), aV, ic: steps.map(s => s.ic), side: steps.map(s => s.side), seg: [seg0, seg1] },
    };
  }

  // ---------- 참고 범위 (초안: 팀 검증 데이터로 확정 예정) ----------
  // 보행 속도: 연령대·성별 평균 ± 0.15 m/s (Bohannon & Andrews 2011 메타분석 값을 근사한 초안, 원문 확인 필요)
  const SPEED_NORM = { M: [1.36, 1.43, 1.43, 1.43, 1.34, 1.26, 0.97], F: [1.34, 1.34, 1.39, 1.31, 1.24, 1.13, 0.94] };
  function refs(subject = {}) {
    const age = subject.age || 30, sex = subject.sex === 'F' ? 'F' : 'M';
    const dec = clamp(Math.floor(age / 10) - 2, 0, 6);
    const sp = SPEED_NORM[sex][dec];
    const decLabel = age >= 80 ? '80대 이상' : `${Math.max(20, Math.floor(age / 10) * 10)}대`;
    return {
      ageGroup: decLabel, sexLabel: sex === 'F' ? '여성' : '남성',
      speed: [+(sp - 0.15).toFixed(2), +(sp + 0.15).toFixed(2)],
      cadence: age >= 70 ? [95, 115] : [100, 120],
      gctBalance: [47.6, 52.4],
      stepLength: subject.height ? [+(subject.height / 100 * 0.40).toFixed(2), +(subject.height / 100 * 0.45).toFixed(2)] : [0.6, 0.75],
      stancePct: [57, 63], doubleSupportPct: [16, 24],
      rotation: [8, 12], obliquity: [6, 10], tilt: [2, 5],
      grfPeak: [1.0, 1.3], si: [0, 10], cv: [0, 4],
    };
  }
  // 등급: 정상(범위 안) / 주의(경계에서 허용폭 이내 이탈) / 심각(허용폭 초과)
  // 허용폭 기본값 = 경계값의 10% (PRD). GCT 밸런스처럼 50% 근처 비율 지표는 절대 허용폭을 따로 둔다.
  const TOL = { gctBalance: 2.4 };
  function grade(v, [lo, hi], key) {
    if (!isFinite(v)) return { level: 'na', label: '-', dev: 0 };
    if (v >= lo && v <= hi) return { level: 'normal', label: '정상', dev: 0 };
    const base = v < lo ? lo : hi;
    const tol = TOL[key] || Math.abs(base) * 0.1 || 1;
    const dev = Math.abs(v - base) / tol * 0.1;   // dev 0.1 = 주의/심각 경계
    const dir = v < lo ? 'low' : 'high';
    return dev <= 0.10 ? { level: 'caution', label: '주의', dev, dir } : { level: 'severe', label: '심각', dev, dir };
  }
  function score(r, ref) {
    const sub = dev => clamp(100 - dev * 300, 0, 100);
    const d = (v, rg, key) => grade(v, rg, key).dev;
    const parts = {
      speed: sub(d(r.general.speed, ref.speed)),
      symmetry: sub(mean([d(r.gct.balance, ref.gctBalance, 'gctBalance'), d(r.symmetry.stepTime, ref.si), d(r.symmetry.stepLength, ref.si)])),
      pelvis: r.pelvis.available ? sub(mean([d(r.pelvis.rotation.rom, ref.rotation), d(r.pelvis.obliquity.rom, ref.obliquity), d(r.pelvis.tilt.rom, ref.tilt)])) : 80,
      rhythm: sub(mean([d(r.general.cadence, ref.cadence), d(r.symmetry.cv, ref.cv)])),
      impact: sub(d(Math.max(r.grf.left, r.grf.right), ref.grfPeak)),
    };
    const total = parts.speed * 0.30 + parts.symmetry * 0.25 + parts.pelvis * 0.20 + parts.rhythm * 0.15 + parts.impact * 0.10;
    return { total: Math.round(total), parts };
  }

  // ---------- 데모용 합성 보행 데이터 ----------
  // 폰 세로, 화면이 뒤(바깥)를 향하게 천골에 고정: 기기 y=위, z=뒤, x=좌우
  function synth(o = {}) {
    const fs = o.fs || 60, rnd = mulberry(o.seed || 7);
    const stepR = o.stepR || 0.52, stepL = o.stepL || 0.58, nSteps = o.steps || 16;
    const pre = 1.4, post = 2.0;
    const ics = []; let tc = pre;
    for (let k = 0; k < nSteps; k++) {
      ics.push({ t: tc, side: k % 2 === 0 ? 'R' : 'L' });
      const base = ics[k].side === 'R' ? stepL : stepR;
      const ramp = k < 2 ? 1.15 - 0.07 * k : k > nSteps - 3 ? 1.08 : 1;
      tc += base * ramp * (1 + (rnd() - 0.5) * 0.04);
    }
    const tEnd = ics[nSteps - 1].t + post;
    const g = 9.80665, samples = [];
    const strideT = stepR + stepL;
    const amp = s => (s === 'L' ? 1.0 : 1.12);
    for (let t = 0; t < tEnd; t += 1 / fs + (rnd() - 0.5) * 0.002) {
      let av = 0, walking = t > pre - 0.3 && t < ics[nSteps - 1].t + 0.5;
      for (const ic of ics) {
        const u = t - ic.t; if (u < -0.3 || u > 0.8) continue;
        const A = amp(ic.side);
        av += A * 3.4 * Math.exp(-0.5 * ((u - 0.10) / 0.045) ** 2) - 2.3 * Math.exp(-0.5 * ((u - 0.30) / 0.085) ** 2) + 0.7 * Math.exp(-0.5 * ((u - 0.47) / 0.05) ** 2);
      }
      av += (rnd() - 0.5) * 0.3;
      const ph = 2 * Math.PI * (t - pre) / strideT, w = walking ? 1 : 0;
      const yawRate = w * 5.0 * (2 * Math.PI / strideT) * Math.cos(ph + 0.9);           // 회전 ±5°
      const oblRate = w * 3.8 * (2 * Math.PI / strideT) * Math.cos(ph + 0.2);           // 경사 ±3.8°
      const tiltRate = w * 1.6 * (4 * Math.PI / strideT) * Math.cos(2 * ph);            // 기울기 ±1.6°
      const ap = w * 1.2 * Math.sin(2 * ph);
      samples.push({
        t: t * 1000,
        ax: (rnd() - 0.5) * 0.2 + w * 0.6 * Math.sin(ph),
        ay: g + av,
        az: -ap + (rnd() - 0.5) * 0.2,
        gx: -tiltRate + (rnd() - 0.5) * 2,
        gy: yawRate + (rnd() - 0.5) * 2,
        gz: -oblRate + (rnd() - 0.5) * 2,
      });
    }
    return { samples, calib: { g: [0, g, 0] }, meta: { synthetic: true } };
  }
  function mulberry(a) { return () => { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }

  function toCSV(raw) {
    const rows = ['t_ms,ax,ay,az,gx,gy,gz'];
    for (const s of raw.samples) rows.push([s.t.toFixed(1), s.ax, s.ay, s.az, s.gx, s.gy, s.gz].map(v => (typeof v === 'number' ? +v.toFixed(4) : v)).join(','));
    return rows.join('\n');
  }
  function fromCSV(text) {
    const lines = text.trim().split(/\r?\n/); const head = lines.shift().split(',');
    const samples = lines.map(l => { const v = l.split(',').map(Number); const o = {}; head.forEach((h, i) => (o[h === 't_ms' ? 't' : h] = v[i])); return o; });
    return { samples };
  }

  const api = { analyze, refs, grade, score, synth, toCSV, fromCSV, version: '0.1.0' };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.GaitCore = api;
})(typeof self !== 'undefined' ? self : this);
