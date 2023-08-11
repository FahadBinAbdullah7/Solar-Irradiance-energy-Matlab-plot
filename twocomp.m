clc
clear all
close all

l = 24.955111; % LATITUDE
n = 1:365; % Days of the year (e.g., 1 to 365)

a = 23.45 * sind((n + 284) * (360 / 365)); % Declination angle

tlt = l; % TILT ANGLE                

AA = 1160 + 75 * sind((360 / 365) * (n - 275));
kk = 0.174 + 0.035 * sind((360 / 365) * (n - 100));
cc = 0.095 + (0.04 * sind((360 / 365) * (n - 100)));
alb = 0.2; % GROUND reflectance

months = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

total_irradiance = zeros(size(n));
total_direct = zeros(size(n));
total_diffusion = zeros(size(n));
total_reflect = zeros(size(n));

monthly_irradiance = zeros(12, 1);
monthly_energy_day15 = zeros(12, 1);

day_counter = 1;

for month = 1:12
    for day = 1:month_days(month)
        Ws = acosd((-tand(l) * tand(a(day_counter)))); % Sunrise angle

        Sr = 12 - ((1 / 15) * (acosd(-tand(l) * tand(a(day_counter))))); % Sunrise time
        Ss = 12 + ((1 / 15) * (acosd(-tand(l) * tand(a(day_counter))))); % Sunset time

        % Rest of the calculations for the current day
        timle = [Sr, floor(Sr) + 1:0.25:floor(Ss), Ss];
        p = length(timle);
        total = zeros(size(timle));
        Idrect = zeros(size(timle));
        idt = zeros(size(timle));
        irt = zeros(size(timle));

        for i = 1:p
            % Hour angle calculation
            ws = (-Ws + (((2 * Ws) / (Ss - Sr)) * (timle(i) - Sr)));
            % Rest of the calculations for the current time of the day
            A = asind((sind(a(day_counter)) * sind(l)) + (cosd(a(day_counter)) * cosd(l) * cosd(ws)));
            Za = 90 - A;
            AM = (1 / cosd(Za));
            AM2 = (1 / sind(A));
            fys = asind((cosd(a(day_counter)) * sind(ws)) / cosd(A));
            kosh = (cosd(A) * cosd(fys - 0) * sind(tlt)) + (sind(A) * cosd(tlt));
            Ib = AA(day_counter) * exp(-kk(day_counter) * AM);

            if (Ib == inf)
                Ib = 0;
            else
                Ib = AA(day_counter) * exp(-kk(day_counter) * AM);
            end

            Idrect(i) = Ib * kosh;
            refactpf = ((1 - cosd(tlt)) / 2);
            difactmf = ((1 + cosd(tlt)) / 2);
            idt(i) = cc(day_counter) * Ib * difactmf;
            irt(i) = alb * Ib * (sind(A) + cc(day_counter)) * refactpf;

            total(i) = irt(i) + idt(i) + Idrect(i);
        end

        total_irradiance(day_counter) = sum(total);
        total_direct(day_counter) = sum(Idrect);
        total_diffusion(day_counter) = sum(idt);
        total_reflect(day_counter) = sum(irt);

        day_counter = day_counter + 1;
    end

    % Calculate the monthly sum of irradiance
    monthly_irradiance(month) = sum(total_irradiance((day_counter - month_days(month)): (day_counter - 1)));
    % Calculate the monthly energy (Day 15 value * total days in month)
    monthly_energy_day15(month) = total_irradiance(day_counter - month_days(month) + 14) * month_days(month);
end

% Define the correct sequence of months
months_ordered = categorical(months, months);

% Plotting both the total energy for each month and the monthly sum of irradiance in the same figure
figure;
subplot(2, 1, 1);
bar(months_ordered, monthly_energy_day15);
xlabel('Months');
ylabel('Total Energy (Watt/m^2)');
title('Total Energy per Month (Using Day 15 Value)');

subplot(2, 1, 2);
bar(months_ordered, monthly_irradiance);
xlabel('Months');
ylabel('Total Irradiance (Watt/m^2)');
title('Total Monthly Irradiance');

% Adjusting the position of the subplot titles
subplot(2, 1, 1);
pos1 = get(gca, 'Position');
pos1(2) = 0.65;
set(gca, 'Position', pos1);

subplot(2, 1, 2);
pos2 = get(gca, 'Position');
pos2(2) = 0.35;
set(gca, 'Position', pos2);

% Adjusting the overall plot title position
sgtitle('Comparison of Total Monthly Irradiance and Total Energy per Month (Using Day 15 Value)');
