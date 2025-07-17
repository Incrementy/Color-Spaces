document.addEventListener('DOMContentLoaded', () => {
    const systemSelect = document.getElementById('color-system-select');
    const randomizeButton = document.getElementById('randomize-button');
    const canvas = document.getElementById('color-canvas');
    const ctx = canvas.getContext('2d');
    const rotationModeSelect = document.getElementById('rotation-mode-select');
    const ditheringSlider = document.getElementById('dithering-slider');
    const ditheringValueDisplay = document.getElementById('dithering-value');

    const currentSystemDisplay = document.getElementById('current-system');
    const currentCenterDisplay = document.getElementById('current-center');
    const currentTiltDisplay = document.getElementById('current-tilt');
    const currentPanDisplay = document.getElementById('current-pan');
    const currentRollDisplay = document.getElementById('current-roll');

    const COLOR_SYSTEMS = {
        rgb: {
            dims: ['R', 'G', 'B'],
            ranges: [[0, 255], [0, 255], [0, 255]],
            toRGB: (c) => c,
            isCylindrical: false
        },
        hsl: {
            dims: ['H', 'S', 'L'],
            ranges: [[0, 360], [0, 100], [0, 100]],
            toRGB: hslToRgb,
            isCylindrical: true
        },
        hsv: {
            dims: ['H', 'S', 'V'],
            ranges: [[0, 360], [0, 100], [0, 100]],
            toRGB: hsvToRgb,
            isCylindrical: true
        },
        lab: {
            dims: ['L*', 'a*', 'b*'],
            ranges: [[0, 100], [-128, 127], [-128, 127]],
            toRGB: labToRgb,
            isCylindrical: false
        },
        lch: {
            dims: ['L*', 'C*', 'h'],
            ranges: [[0, 100], [0, 134], [0, 360]],
            toRGB: lchToRgb,
            isCylindrical: true
        },
        cmyk: {
            dims: ['C', 'M', 'Y', 'K'],
            ranges: [[0, 100], [0, 100], [0, 100]],
            toRGB: cmykToRgb,
            isCylindrical: false
        },
        xyz: {
            dims: ['X', 'Y', 'Z'],
            ranges: [[0, 95.047], [0, 100], [0, 108.883]],
            toRGB: xyzToRgb,
            isCylindrical: false
        },
        yiq: {
            dims: ['Y', 'I', 'Q'],
            ranges: [[0, 1], [-0.5957, 0.5957], [-0.5226, 0.5226]],
            toRGB: yiqToRgb,
            isCylindrical: false
        },
        oklab: {
            dims: ['L', 'a', 'b'],
            ranges: [[0, 1], [-0.4, 0.4], [-0.4, 0.4]],
            toRGB: oklabToRgb,
            isCylindrical: false
        },
    };

    let currentState = {
        system: 'rgb',
        centerPoint: [0, 0, 0],
        tilt: 0,
        pan: 0,
        roll: 0,
        rotationMode: 'none',
        ditheringLevel: 0
    };

    let lastTime = 0;
    const rotationSpeedRoll = 30;
    const rotationSpeedPitch = 20;
    const rotationSpeedYaw = 25;

    const vec3 = {
        create: (x = 0, y = 0, z = 0) => [x, y, z],
        add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
        scale: (v, s) => [v[0] * s, v[1] * s, v[2] * s],
        dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
        cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
        normalize: (v) => {
            const len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return len > 1e-6 ? [v[0] / len, v[1] / len, v[2] / len] : [0, 0, 0];
        }
    };

    const mat3 = {
        identity: () => [1, 0, 0, 0, 1, 0, 0, 0, 1],
        rotateX: (rad) => {
            const c = Math.cos(rad);
            const s = Math.sin(rad);
            return [1, 0, 0, 0, c, s, 0, -s, c];
        },
        rotateY: (rad) => {
            const c = Math.cos(rad);
            const s = Math.sin(rad);
            return [c, 0, -s, 0, 1, 0, s, 0, c];
        },
        rotateZ: (rad) => {
            const c = Math.cos(rad);
            const s = Math.sin(rad);
            return [c, s, 0, -s, c, 0, 0, 0, 1];
        },
        multiply: (a, b) => {
            let result = [];
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    result[i * 3 + j] = 0;
                    for (let k = 0; k < 3; k++) {
                        result[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
                    }
                }
            }
            return result;
        },
        multiplyVec3: (m, v) => [
            m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
            m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
            m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
        ]
    };

    function generateRandomParameters() {
        const system = COLOR_SYSTEMS[currentState.system];
        const reductionFactor = 0.9; 

        currentState.centerPoint = system.ranges.map((range, index) => {
            let min = range[0];
            let max = range[1];
            let fullRange = max - min;
            let reducedRange = fullRange * reductionFactor;
            let offset = (fullRange - reducedRange) / 2;

            let newMin = min + offset;
            let newMax = max - offset;
            
            if (system.isCylindrical && index === 0) {
                return range[0] + Math.random() * (range[1] - range[0]);
            } else {
                return newMin + Math.random() * (newMax - newMin);
            }
        });

        if (currentState.system === 'cmyk') {
            const kRange = system.ranges[3];
            const kFullRange = kRange[1] - kRange[0];
            const kReducedRange = kFullRange * reductionFactor;
            const kOffset = (kFullRange - kReducedRange) / 2;
            currentState.centerPoint[3] = (kRange[0] + kOffset) + Math.random() * kReducedRange;
        } else if (currentState.system === 'yiq') {
            currentState.centerPoint[0] = 0.5 + (Math.random() * 0.4 - 0.2);
        } else if (currentState.system === 'oklab') {
            currentState.centerPoint[0] = 0.5 + (Math.random() * 0.4 - 0.2);
        }

        currentState.tilt = Math.random() * 180 - 90;
        currentState.pan = Math.random() * 360;
        currentState.roll = Math.random() * 360;

        currentState.initialTilt = currentState.tilt;
        currentState.initialPan = currentState.pan;
        currentState.initialRoll = currentState.roll;

        updateDisplayedParameters();
    }

    function updateDisplayedParameters() {
        const system = COLOR_SYSTEMS[currentState.system];
        currentSystemDisplay.textContent = system.dims.slice(0, 3).join('/');
        let centerText = currentState.centerPoint.slice(0, 3).map(val => val.toFixed(3)).join(', ');
        if (currentState.system === 'cmyk') {
            centerText = currentState.centerPoint.map((val, i) => i === 3 ? `K:${val.toFixed(1)}` : val.toFixed(1)).join(', ');
        }
        currentCenterDisplay.textContent = centerText;
        currentTiltDisplay.textContent = `${currentState.tilt.toFixed(1)}°`;
        currentPanDisplay.textContent = `${currentState.pan.toFixed(1)}°`;
        currentRollDisplay.textContent = `${currentState.roll.toFixed(1)}°`;
    }

    function applyDithering(data, width, height, ditherAmount) {
        if (ditherAmount === 0) return;

        const maxColor = 255;
        const quantize = (val, max) => Math.round(val / max) * max;

        const errors = new Float32Array(width * height * 3);

        for (let y = 0; y < height; y++) {
            for (let x = 0; x < width; x++) {
                const i = (y * width + x) * 4;
                const errorIdx = (y * width + x) * 3;

                let oldR = data[i] + errors[errorIdx];
                let oldG = data[i + 1] + errors[errorIdx + 1];
                let oldB = data[i + 2] + errors[errorIdx + 2];

                let newR = quantize(oldR, maxColor);
                let newG = quantize(oldG, maxColor);
                let newB = quantize(oldB, maxColor);
                
                data[i] = newR;
                data[i + 1] = newG;
                data[i + 2] = newB;

                let errR = (oldR - newR) * (ditherAmount / 100);
                let errG = (oldG - newG) * (ditherAmount / 100);
                let errB = (oldB - newB) * (ditherAmount / 100);

                if (x + 1 < width) {
                    errors[((y * width) + (x + 1)) * 3] += errR * (7 / 16);
                    errors[((y * width) + (x + 1)) * 3 + 1] += errG * (7 / 16);
                    errors[((y * width) + (x + 1)) * 3 + 2] += errB * (7 / 16);
                }
                if (x > 0 && y + 1 < height) {
                    errors[(((y + 1) * width) + (x - 1)) * 3] += errR * (3 / 16);
                    errors[(((y + 1) * width) + (x - 1)) * 3 + 1] += errG * (3 / 16);
                    errors[(((y + 1) * width) + (x - 1)) * 3 + 2] += errB * (3 / 16);
                }
                if (y + 1 < height) {
                    errors[(((y + 1) * width) + x) * 3] += errR * (5 / 16);
                    errors[(((y + 1) * width) + x) * 3 + 1] += errG * (5 / 16);
                    errors[(((y + 1) * width) + x) * 3 + 2] += errB * (5 / 16);
                }
                if (x + 1 < width && y + 1 < height) {
                    errors[(((y + 1) * width) + (x + 1)) * 3] += errR * (1 / 16);
                    errors[(((y + 1) * width) + (x + 1)) * 3 + 1] += errG * (1 / 16);
                    errors[(((y + 1) * width) + (x + 1)) * 3 + 2] += errB * (1 / 16);
                }
            }
        }
    }


    function drawArbitraryCrossSection() {
        const { width, height } = canvas;
        const system = COLOR_SYSTEMS[currentState.system];
        const imageData = ctx.createImageData(width, height);
        const data = imageData.data;

        const center = currentState.centerPoint;
        const tiltRad = currentState.tilt * Math.PI / 180;
        const panRad = currentState.pan * Math.PI / 180;
        const rollRad = currentState.roll * Math.PI / 180;

        let rotationMatrix = mat3.identity();
        rotationMatrix = mat3.multiply(mat3.rotateY(panRad), rotationMatrix);
        rotationMatrix = mat3.multiply(mat3.rotateX(tiltRad), rotationMatrix);
        rotationMatrix = mat3.multiply(mat3.rotateZ(rollRad), rotationMatrix);

        const initialPlaneXAxis = vec3.create(1, 0, 0);
        const initialPlaneYAxis = vec3.create(0, 1, 0);

        const planeXAxis = vec3.normalize(mat3.multiplyVec3(rotationMatrix, initialPlaneXAxis));
        const planeYAxis = vec3.normalize(mat3.multiplyVec3(rotationMatrix, initialPlaneYAxis));

        let scaleFactor;
        if (currentState.system === 'hsl' || currentState.system === 'hsv') {
            scaleFactor = 60;
        } else if (currentState.system === 'lab' || currentState.system === 'lch') {
            scaleFactor = 110;
        } else if (currentState.system === 'oklab') {
            scaleFactor = 0.25;
        } else {
            const maxRange = system.ranges.slice(0, 3).reduce((max, range) => Math.max(max, range[1] - range[0]), 0);
            scaleFactor = maxRange / 2;
        }

        const pixelDataForDithering = new Float32Array(width * height * 4);

        for (let y = 0; y < height; y++) {
            for (let x = 0; x < width; x++) {
                const u = (x / (width - 1)) * 2 - 1;
                const v = (1 - y / (height - 1)) * 2 - 1;

                const point3D = vec3.add(
                    center.slice(0, 3),
                    vec3.add(
                        vec3.scale(planeXAxis, u * scaleFactor),
                        vec3.scale(planeYAxis, v * scaleFactor)
                    )
                );

                let colorSpacePoint = [...point3D];

                if (currentState.system === 'cmyk') {
                    colorSpacePoint.push(currentState.centerPoint[3]); 
                }
                
                const clampedPoint = colorSpacePoint.map((val, i) => {
                    const [min, max] = system.ranges[i] || [val, val];
                    return Math.max(min, Math.min(max, val));
                });

                const [r, g, b] = system.toRGB(clampedPoint);

                const index = (y * width + x) * 4;
                pixelDataForDithering[index] = r;
                pixelDataForDithering[index + 1] = g;
                pixelDataForDithering[index + 2] = b;
                pixelDataForDithering[index + 3] = 255;
            }
        }

        for(let i = 0; i < width * height * 4; i++) {
            data[i] = pixelDataForDithering[i];
        }

        applyDithering(data, width, height, currentState.ditheringLevel);

        ctx.putImageData(imageData, 0, 0);
    }

    function animate(currentTime) {
        if (!lastTime) lastTime = currentTime;
        const deltaTime = (currentTime - lastTime) / 1000;
        lastTime = currentTime;

        let shouldRedraw = false;

        switch (currentState.rotationMode) {
            case 'pitch':
                currentState.tilt = (currentState.tilt + rotationSpeedPitch * deltaTime) % 360;
                shouldRedraw = true;
                break;
            case 'yaw':
                currentState.pan = (currentState.pan + rotationSpeedYaw * deltaTime) % 360;
                shouldRedraw = true;
                break;
            case 'roll':
                currentState.roll = (currentState.roll + rotationSpeedRoll * deltaTime) % 360;
                shouldRedraw = true;
                break;
            case 'pitch_yaw':
                currentState.tilt = (currentState.tilt + rotationSpeedPitch * deltaTime * 0.8) % 360;
                currentState.pan = (currentState.pan + rotationSpeedYaw * deltaTime * 0.8) % 360;
                shouldRedraw = true;
                break;
            case 'pitch_roll':
                currentState.tilt = (currentState.tilt + rotationSpeedPitch * deltaTime * 0.8) % 360;
                currentState.roll = (currentState.roll + rotationSpeedRoll * deltaTime * 0.8) % 360;
                shouldRedraw = true;
                break;
            case 'yaw_roll':
                currentState.pan = (currentState.pan + rotationSpeedYaw * deltaTime * 0.8) % 360;
                currentState.roll = (currentState.roll + rotationSpeedRoll * deltaTime * 0.8) % 360;
                shouldRedraw = true;
                break;
            case 'spiral':
                const spiralSpeedMultiplier = 1.5;
                const spiralTimeFactor = (currentTime / 1000) * spiralSpeedMultiplier; 

                currentState.pan = (currentState.initialPan + spiralTimeFactor * rotationSpeedYaw * 0.5) % 360;
                currentState.tilt = (currentState.initialTilt + Math.sin(spiralTimeFactor * 0.2) * 45); 
                currentState.tilt = Math.max(-90, Math.min(90, currentState.tilt));
                currentState.roll = (currentState.initialRoll + spiralTimeFactor * rotationSpeedRoll * 0.05) % 360;

                shouldRedraw = true;
                break;
            case 'none':
            default:
                break;
        }
        
        if (shouldRedraw) {
            drawArbitraryCrossSection();
            updateDisplayedParameters();
        }

        requestAnimationFrame(animate);
    }

    systemSelect.addEventListener('change', (e) => {
        currentState.system = e.target.value;
        generateRandomParameters();
        drawArbitraryCrossSection();
        lastTime = 0; 
    });

    randomizeButton.addEventListener('click', () => {
        generateRandomParameters();
        drawArbitraryCrossSection();
        lastTime = 0; 
    });

    rotationModeSelect.addEventListener('change', (e) => {
        currentState.rotationMode = e.target.value;
        lastTime = 0; 
        currentState.initialTilt = currentState.tilt;
        currentState.initialPan = currentState.pan;
        currentState.initialRoll = currentState.roll;
        drawArbitraryCrossSection();
    });

    ditheringSlider.addEventListener('input', (e) => {
        currentState.ditheringLevel = parseInt(e.target.value, 10);
        updateDitheringDisplay(currentState.ditheringLevel);
        drawArbitraryCrossSection();
    });

    function clampRgb(linearRgbArray) {
        let r_linear = linearRgbArray[0];
        let g_linear = linearRgbArray[1];
        let b_linear = linearRgbArray[2];

        return [r_linear, g_linear, b_linear].map(v => Math.max(0, Math.min(1, linearRgbToSrgb(v))) * 255);
    }

    function hslToRgb([h, s, l]) {
        s /= 100; l /= 100;
        let c = (1 - Math.abs(2 * l - 1)) * s,
            x = c * (1 - Math.abs((h / 60) % 2 - 1)),
            m = l - c / 2,
            r_linear = 0, g_linear = 0, b_linear = 0; 

        if (0 <= h && h < 60) { [r_linear, g_linear, b_linear] = [c, x, 0]; }
        else if (60 <= h && h < 120) { [r_linear, g_linear, b_linear] = [x, c, 0]; }
        else if (120 <= h && h < 180) { [r_linear, g_linear, b_linear] = [0, c, x]; }
        else if (180 <= h && h < 240) { [r_linear, g_linear, b_linear] = [0, x, c]; }
        else if (240 <= h && h < 300) { [r_linear, g_linear, b_linear] = [x, 0, c]; }
        else if (300 <= h && h < 360) { [r_linear, g_linear, b_linear] = [c, 0, x]; }
        
        return clampRgb([r_linear + m, g_linear + m, b_linear + m]);
    }

    function hsvToRgb([h, s, v]) {
        s /= 100; v /= 100;
        let c = v * s,
            x = c * (1 - Math.abs((h / 60) % 2 - 1)),
            m = v - c,
            r_linear = 0, g_linear = 0, b_linear = 0; 

        if (0 <= h && h < 60) { [r_linear, g_linear, b_linear] = [c, x, 0]; }
        else if (60 <= h && h < 120) { [r_linear, g_linear, b_linear] = [x, c, 0]; }
        else if (120 <= h && h < 180) { [r_linear, g_linear, b_linear] = [0, c, x]; }
        else if (180 <= h && h < 240) { [r_linear, g_linear, b_linear] = [0, x, c]; }
        else if (240 <= h && h < 300) { [r_linear, g_linear, b_linear] = [x, 0, c]; }
        else if (300 <= h && h < 360) { [r_linear, g_linear, b_linear] = [c, 0, x]; }

        return clampRgb([r_linear + m, g_linear + m, b_linear + m]);
    }

    function lchToRgb([l, c, h]) {
        const hRad = h * (Math.PI / 180);
        const a = c * Math.cos(hRad);
        const b = c * Math.sin(hRad);
        return labToRgb([l, a, b]);
    }

    function labToRgb([l, a, b]) {
        const Xn = 95.047;
        const Yn = 100.000;
        const Zn = 108.883;

        let fy = (l + 16) / 116;
        let fx = a / 500 + fy;
        let fz = fy - b / 200;

        const delta = 6 / 29;
        let x_xyz = fx > delta ? fx * fx * fx : (fx - 16 / 116) * 3 * delta * delta;
        let y_xyz = fy > delta ? fy * fy * fy : (fy - 16 / 116) * 3 * delta * delta;
        let z_xyz = fz > delta ? fz * fz * fz : (fz - 16 / 116) * 3 * delta * delta;

        x_xyz *= Xn;
        y_xyz *= Yn;
        z_xyz *= Zn;

        x_xyz /= 100; y_xyz /= 100; z_xyz /= 100;

        let r_linear = x_xyz * 3.2406 + y_xyz * -1.5372 + z_xyz * -0.4986;
        let g_linear = x_xyz * -0.9689 + y_xyz * 1.8758 + z_xyz * 0.0415;
        let b_linear = x_xyz * 0.0557 + y_xyz * -0.2040 + z_xyz * 1.0570;

        return clampRgb([r_linear, g_linear, b_linear]);
    }

    function cmykToRgb([c, m, y, k]) {
        c /= 100;
        m /= 100;
        y /= 100;
        k /= 100;

        let r_linear = 1 - Math.min(1, c * (1 - k) + k);
        let g_linear = 1 - Math.min(1, m * (1 - k) + k);
        let b_linear = 1 - Math.min(1, y * (1 - k) + k);

        return clampRgb([r_linear, g_linear, b_linear]);
    }

    function xyzToRgb([x, y, z]) {
        x /= 100; y /= 100; z /= 100;

        let r_linear = x * 3.2406 + y * -1.5372 + z * -0.4986;
        let g_linear = x * -0.9689 + y * 1.8758 + z * 0.0415;
        let b_linear = x * 0.0557 + y * -0.2040 + z * 1.0570;

        return clampRgb([r_linear, g_linear, b_linear]);
    }

    function yiqToRgb([y, i, q]) {
        let r_linear = y + 0.956 * i + 0.621 * q;
        let g_linear = y - 0.272 * i - 0.647 * q;
        let b_linear = y - 1.106 * i + 1.703 * q;
        
        return clampRgb([r_linear, g_linear, b_linear]);
    }

    const D65 = [0.95047, 1.00000, 1.08883]; 

    function srgbToLinearRgb(c) {
        return c > 0.04045 ? Math.pow((c + 0.055) / 1.055, 2.4) : c / 12.92;
    }

    function linearRgbToSrgb(c) {
        return c > 0.0031308 ? (1.055 * Math.pow(c, 1 / 2.4) - 0.055) : 12.92 * c;
    }

    function linearRgbToXyz([r, g, b]) {
        r = srgbToLinearRgb(r / 255);
        g = srgbToLinearRgb(g / 255);
        b = srgbToLinearRgb(b / 255);

        const x = r * 0.4124564 + g * 0.3575761 + b * 0.1804375;
        const y = r * 0.2126729 + g * 0.7151522 + b * 0.0721750;
        const z = r * 0.0193339 + g * 0.1191920 + b * 0.9503041;
        return [x * 100, y * 100, z * 100];
    }

    function xyzToLinearRgb([x, y, z]) {
        x /= 100; y /= 100; z /= 100;

        let r = x * 3.2406 + y * -1.5372 + z * -0.4986;
        let g = x * -0.9689 + y * 1.8758 + z * 0.0415;
        let b_linear = x * 0.0557 + y * -0.2040 + z * 1.0570;

        return clampRgb([r, g, b_linear]);
    }

    const CUBE_ROOT = (x) => Math.cbrt(x);
    const POW3 = (x) => x * x * x;

    function oklabToRgb([l, a, b]) {
        const l_ = l + 0.3963377774 * a + 0.2158037573 * b;
        const m_ = l - 0.1055613423 * a - 0.0638541728 * b;
        const s_ = l - 0.0894841775 * a - 1.2914855480 * b;

        const L = POW3(l_);
        const M = POW3(m_);
        const S = POW3(s_);

        const r_linear = +4.0767416621 * L - 3.3077115913 * M + 0.2309699292 * S;
        const g_linear = -1.2684380046 * L + 2.6097574011 * M - 0.3413193965 * S;
        const b_linear = -0.0041960863 * L - 0.7034186147 * M + 1.7076147010 * S;
        
        return clampRgb([r_linear, g_linear, b_linear]);
    }

    function init() {
        currentState.system = systemSelect.value;
        currentState.ditheringLevel = parseInt(ditheringSlider.value, 10);
        updateDitheringDisplay(currentState.ditheringLevel);
        generateRandomParameters();
        drawArbitraryCrossSection();
        requestAnimationFrame(animate);
    }

    function updateDitheringDisplay(value) {
        if (value === 0) {
            ditheringValueDisplay.textContent = 'None';
        } else if (value === 100) {
            ditheringValueDisplay.textContent = 'High';
        } else {
            ditheringValueDisplay.textContent = value + '%';
        }
    }

    init();
});