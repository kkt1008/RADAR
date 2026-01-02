function radarSNR_GUI
    % Range vs Radar SNR

    % figure 생성
    % Name 속성(제목) : Radar SNR Calculator
    % GUI window의 위치와 크기 : (100, 100) & 700x450
    fig = uifigure('Name', 'Radar SNR Calculator', 'Position', [100, 100, 700, 450]);

    % UI 요소들 나열(입력 Field 옆에 표시될 예정)
    labels = ["Target Min Range (m)", "Target Max Range (m)", ...
              "System Temperature (K)", "Bandwidth (Hz)", ...
              "Noise Figure (dB)", "Losses (dB)", ...
              "Peak Power (W)", "Antenna Gain (dB)", ...
              "Frequency (Hz)", "Target RCS (dBsm)"];
          
    % 수치 입력 필드 배열 초기화
    fields = gobjects(1, 10);

    % for 반복문을 활용해 Field 생성
    for i = 1:10
        uilabel(fig, 'Text', labels(i), 'Position', [20, 400 - 30*i, 150, 22]);
        fields(i) = uieditfield(fig, 'numeric', ...
                                'Position', [180, 400 - 30*i, 100, 22]);
    end

    % 그래프 plot 좌표 설명 및 세부 정보 표시
    ax = uiaxes(fig, 'Position', [300, 60, 380, 360]);
    title(ax, 'Output SNR vs Target Range');
    xlabel(ax, 'Target Range (km)');
    ylabel(ax, 'Output SNR (dB)');
    grid(ax, 'on');

    % Button 생성
    btn = uibutton(fig, 'Text', 'Calculate SNR', ...
                   'Position', [80, 60, 160, 30], ...
                   'ButtonPushedFcn', @(btn,event) plotSNR());

    % Callback function
    function plotSNR()
        try
            % Get input values from GUI
            minR = fields(1).Value;
            maxR = fields(2).Value;
            T = fields(3).Value;
            B = fields(4).Value;
            NF_dB = fields(5).Value;
            L_dB = fields(6).Value;
            Pt = fields(7).Value;
            G_dB = fields(8).Value;
            f = fields(9).Value;
            rcs_dBsm = fields(10).Value;

            % Validation
            if any(isnan([minR, maxR, T, B, NF_dB, L_dB, Pt, G_dB, f, rcs_dBsm]))
                uialert(fig, 'All input fields must be filled.', 'Missing Input');
                return;
            end

            % Constants
            k = 1.38e-23;
            c = 3e8;
            lambda = c / f;

            % Convert dB values to linear
            NF = 10^(NF_dB / 10);
            L = 10^(L_dB / 10);
            G = 10^(G_dB / 10);
            sigma = 10^(rcs_dBsm / 10);

            % Compute noise power
            N = k * T * B * NF;

            % Range Array 생성 (최소거리 ~ 최대거리 with element 1000개)
            R = linspace(minR, maxR, 1000);

            % Radar equation
            Pr = (Pt * G^2 * lambda^2 * sigma) ./ ((4 * pi)^3 * R.^4 * L);
            SNR = Pr / N;
            SNR_dB = 10 * log10(SNR);

            % Plot
            plot(ax, R / 1000, SNR_dB, 'k', 'LineWidth', 1);
            title(ax, 'Output Signal to Noise Ratio vs Target Range');
            xlabel(ax, 'Target Range (km)');
            ylabel(ax, 'Output SNR (dB)');
            grid(ax, 'on');
        catch ME
            uialert(fig, ME.message, 'Error');
        end
    end
end
