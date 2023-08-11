clc
clear all
close all

l = 24.95511; % LATITUDE
n = 1:365; % Days of the year ( 1 to 365)

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

average_daily_energy = zeros(12, 1);

day_counter = 1;

for month = 1:12
    total_monthly_irradiance = 0; % Initialize the sum of irradiance for the current month
    
    for day = 1:month_days(month)
        Ws = acosd((-tand(l) * tand(a(day_counter)))); % Sunrise angle

        Sr = 12 - ((1 / 15) * (acosd(-tand(l) * tand(a(day_counter))))); % Sunrise time
        Ss = 12 + ((1 / 15) * (acosd(-tand(l) * tand(a(day_counter))))); % Sunset time

        %  calculations for the current day
        timle = [Sr, floor(Sr) + 1:0.25:floor(Ss), Ss];
        p = length(timle);
        total = zeros(size(timle));
        Idrect = zeros(size(timle));
        idt = zeros(size(timle));
        irt = zeros(size(timle));

        for i = 1:p
            % Hour angle calculation
            ws = (-Ws + (((2 * Ws) / (Ss - Sr)) * (timle(i) - Sr)));
            %  current time of the day
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
        
        % Adding daily irradiance to the sum for the current month
        total_monthly_irradiance = total_monthly_irradiance + sum(total);
        
        day_counter = day_counter + 1;
    end

    % Calculate the average daily energy for the current month
    average_daily_energy(month) = total_monthly_irradiance / month_days(month);
end

%  sequence of months in the categorical array
months_ordered = categorical(months, months);

% Plotting the average daily energy for each month as a bar chart
figure;
bar(months_ordered, average_daily_energy);
xlabel('Months');
ylabel('Average Daily Energy (Watt/m^2)');
title('Average Daily Energy per Month');
