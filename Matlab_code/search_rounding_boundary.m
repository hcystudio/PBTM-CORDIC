%% Search Rounding Boundary
%  21-Sep-2025 14:24:10 @ hcy
%  

clc;
clear;
close all;

%% User Config
k = 9; % Mapping Depth
iter_max = 28; % M
w_frac = 41;

%% Prepare
dd_numax = 2^k-1;
atan_table = zeros(1, 200);
Bn_table = zeros(1, 200);
for i = 1:200
    atan_table(i) = atan(2^(-(i-1)));
end
Bn_table(1) = 1;
for i = 2:200
    Bn_table(i) =  prod((1+2.^(-2*(1:i-1))).*cos(atan_table(2:i))) ;
end

%% Fixed Design
w_phi = iter_max+1;
p_upper = fi([], 1, w_frac+2, w_frac);
p_comp_lower = fi([], 0, w_frac+1, w_frac);
p_Bn_table = fi([], 0, w_frac+1, w_frac);
p_atan_table = fi([], 0, w_phi, w_phi);
p_current_phi = fi([], 1, w_phi+2, w_phi);
cast_Bn_table = cast(Bn_table, 'like', p_Bn_table);
cast_atan_table = cast(atan_table, 'like', p_atan_table);
initial_phi = cast_atan_table(1);
ppp = {p_upper p_comp_lower p_current_phi};

%% Binary Traversal
ddw = ceil(log2(dd_numax)) + 1;
range_up = zeros(1, dd_numax+1);
range_down = zeros(1, dd_numax+1);
ddidx_array = zeros(1, dd_numax+1);
for di = 0:dd_numax
ddidx = ['0' dec2bin(di, ddw - 1)];

ddd = zeros(1, ddw);
for ii = 1:ddw
    if(ddidx(ii) == '1')
        ddd(ii) = -1;
    else
        ddd(ii) = 1;
    end
end

acc_phi = sum(atan_table(1:ddw).*ddd(1:ddw));
last_phi = sum(atan_table(ddw+1:iter_max));
range_up(di+1) = acc_phi + last_phi;
range_down(di+1) = acc_phi - last_phi;

% valid bond
start_up = cos(range_up(di+1));
target_lower = cast(start_up, 'like', p_upper);
[current_phi, dd] = invc_acos(target_lower, iter_max, ppp, cast_atan_table, cast_Bn_table);

if(strcmp(dd(1:k+1), ddidx) == 0)
    fprintf(['Target dd(' ddidx ') upper bound not match! use ' dd(1:k+1) '\n']);
end

start_up = cos(range_down(di+1));
target_lower = cast(start_up, 'like', p_upper);
[current_phi, dd] = invc_acos(target_lower, iter_max, ppp, cast_atan_table, cast_Bn_table);

if(strcmp(dd(1:k+1), ddidx) == 0)
    fprintf(['Target dd(' ddidx ') lower bound not match! use ' dd(1:k+1) '\n']);
end

end

save(['ranges_k' num2str(k) '.mat'], "range_up", "range_down");
% load(['ranges_k' num2str(k) '.mat'], "range_up", "range_down");

figure;
hold on;
scatter(0:dd_numax, range_up, ones(1, length(range_up))*8, 'filled');
scatter(0:dd_numax, range_down, ones(1, length(range_up))*8, 'filled');
figure;
hold on;
scatter(0:dd_numax, cos(range_up), ones(1, length(range_up))*8, 'filled');
scatter(0:dd_numax, cos(range_down), ones(1, length(range_up))*8, 'filled');


range_up_fcos = cast(cos(range_up), 'like', p_upper);
range_dn_fcos = cast(cos(range_down), 'like', p_upper);
%% Make Rounding Boundary Table
intlv_table = zeros(1, 1000);
intlv_table_fixed = cast(zeros(1, 1000), 'like', p_upper);
dd_table = zeros(1, 1000);


start_num = sum(cos(range_up) <= 0)-1;
[~, ui] = max(cos(range_down));
end_num = ui-1;
intlv_table(1) = 0;
intlv_table_fixed(1) = 0;
str_mid = find_rounding_boundary(range_dn_fcos(start_num+1), range_up_fcos(start_num+2));
num_mid = bin2dec(str_mid)/2^p_upper.FractionLength;

intlv_table(2) = num_mid; % start_num
intlv_table_fixed(2) = cast(num_mid, 'like', p_upper);
dd_table(1) = start_num;
tcount = 3;
di = start_num+1;
% debug_ud_un_count = 0;
while di <= end_num-1
    temp_i = di+2;
    while(range_dn_fcos(temp_i) < range_dn_fcos(di+1))
        temp_i = temp_i + 1;
    end

    str_mid = find_rounding_boundary(range_dn_fcos(di+1), range_up_fcos(temp_i));
    num_mid = bin2dec(str_mid)/2^p_upper.FractionLength;

    intlv_table(tcount) = num_mid;
    intlv_table_fixed(tcount) = cast(num_mid, 'like', p_upper);
    dd_table(tcount-1) = di;
    tcount = tcount + 1;
    di = temp_i - 1;
end

intlv_table(tcount) = 1;
intlv_table_fixed(tcount) = 1;
dd_table(tcount-1) = end_num;

intlv_table(tcount+1:end) = '';
intlv_table_fixed(tcount+1:end) = '';
dd_table(tcount:end) = '';


figure;
plot(intlv_table);
figure;
plot(dd_table);
dd_table_new = dd_table;
intlv_table_new = intlv_table;
save(['invc_param_k' num2str(k) '.mat'], "dd_table", "intlv_table");
% load("invc_tab_param.mat", "dd_table", "intlv_table");

for ii = 1:length(intlv_table_fixed)
    fprintf([bin(intlv_table_fixed(ii)) ' ' num2str(ii-1)  '\n']);
end



function [str_mid] = find_rounding_boundary(fnum1, fnum2)
    str1 = bin(fnum1);
    str2 = bin(fnum2);
    
    for i  = 1:length(str1)
        if(strcmp(str1(i), str2(i)))
            continue;
        else
            if(i == length(str1))
                str_mid = str1;
                return;
            end
            str_mid = [str2(1:i-1) '10' repmat('0', [1 length(str2)-i-1])];
            return;
        end
    end
    str_mid = str1;
end


function [current_phi, dd] = invc_acos(target_lower, iter, ppp, cast_atan_table, cast_Bn_table)

    p_upper = ppp{1};
    p_comp_lower = ppp{2};
    p_current_phi = ppp{3};
    % current_lower = 1;
    lower = cast(2^(-1/2), 'like', p_upper);
    upper = cast(2^(-1/2), 'like', p_upper);
    current_phi = cast_atan_table(1);
    dd(1) = '0';
    for i = 2:iter
        tu = upper;
        tl = lower;
        atan_tab = cast_atan_table(i);
        mi = i - 1;
        if(i < 17)
            comp_lower = cast(target_lower * cast_Bn_table(i-1), 'like', p_comp_lower);
        else
            comp_lower = cast(target_lower * cast_Bn_table(end), 'like', p_comp_lower);
        end
    
        la = lower >= comp_lower;
        lb = lower < comp_lower;
        tab_inv = 1;
        if(tu < 0)
            la = ~ la;
            lb = ~ lb;
            tab_inv = -1;
        end
        if(la)
            % inv clkwise
            upper = cast(tu + bitshift(tl, -mi), 'like', p_upper);
            lower = cast(tl - bitshift(tu, -mi), 'like', p_upper);
    
            current_phi = cast(current_phi + atan_tab, 'like', p_current_phi);
            dd(i) = '0';
        elseif(lb)
            % clkwise
            upper = cast(tu - bitshift(tl, -mi), 'like', p_upper);
            lower = cast(tl + bitshift(tu, -mi), 'like', p_upper);

            current_phi = cast(current_phi - atan_tab, 'like', p_current_phi);
            dd(i) = '1';
        else
            fprintf('break at iter:%d\n', i);
            break;
        end
    end
%     upper = abs(upper);
    current_phi = abs(current_phi);
end