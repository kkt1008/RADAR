function maxDetectionRange_GUI
    % 최대 탐지 거리 계산 GUI 만듦

    fig = uifigure('Name', 'Max Detection Range Calculator', ...
        'Position', [100, 100, 500, 420]);

    % 입력 라벨 및 필드 생성
    labels = ["Peak Power (W)", "Antenna Gain (dB)", "Frequency (Hz)", ...
              "Target RCS (dBsm)", "System Temperature (K)", "Bandwidth (Hz)", ...
              "Noise Figure (dB)", "Losses (dB)", "Required SNR (dB)"];
    fields = gobjects(1, numel(labels));

    for i = 1:numel(labels)
        % 각 입력 필드 위치 설정
        uilabel(fig, 'Text', labels(i), 'Position', [20, 400 - 35*i, 180, 22]);
        fields(i) = uieditfield(fig, 'numeric', ...
                                'Position', [210, 400 - 35*i, 120, 22], ...
                                'Placeholder', 'Enter vZalue...');
    end

    % 결과 표시용 라벨
    resultLabel = uilabel(fig, 'Text', 'Max Detection Range: ', ...
                          'Position', [20, 50, 460, 22], ...z
                          'FontWeight', 'bold');

    % 계산 버튼 생성
    uibutton(fig, 'Text', 'Calculate Max Range', ...
        'Position', [150, 15, 200, 30], ...
        'ButtonPushedFcn', @(btn, event) calculateRange());

    % 계산 함수
    function calculateRange()
        try
            % 입력값 가져옴
            Pt = fields(1).Value;
            G_dB = fields(2).Value;
            f = fields(3).Value;
            rcs_dBsm = fields(4).Value;
            T = fields(5).Value;
            B = fields(6).Value;
            NF_dB = fields(7).Value;
            L_dB = fields(8).Value;
            SNR_dB = fields(9).Value;

            % 입력 안 한 필드 있으면 경고
            if any(isnan([Pt, G_dB, f, rcs_dBsm, T, B, NF_dB, L_dB, SNR_dB]))
                uialert(fig, 'Please fill in all input fields.', 'Missing Input');
                return;
            end

            % 상수 및 계산용 변수 설정
            k = 1.38e-23;    % 볼츠만 상수
            c = 3e8;         % 빛 속도
            lambda = c / f;  % 파장

            % dB -> 선형 변환
            G = 10^(G_dB / 10);
            sigma = 10^(rcs_dBsm / 10);
            NF = 10^(NF_dB / 10);
            L = 10^(L_dB / 10);
            SNR = 10^(SNR_dB / 10);

            % 노이즈 전력 계산
            N = k * T * B * NF;

            % 최대 탐지 거리 계산식 (레이더 방정식)
            R = ((Pt * G^2 * lambda^2 * sigma) / ((4*pi)^3 * N * SNR * L))^(1/4);

            % 결과 표시
            resultLabel.Text = sprintf('Max Detection Range: %.2f meters (%.2f km)', R, R/1000);
        catch ME
            uialert(fig, ME.message, 'Error');
        end
    end
end
