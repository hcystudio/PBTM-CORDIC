%% Construct PBT-k of CORDIC
%  2025-09-23 16:18:08 @ hcy
%  Inverse Rotation Mode, circular system

clc;
clear;
close all;

%% User Config
M = 32; % Iteration Number
k = 9;  % Mapping Depth
w_frac = 41; % Fraction Width
arctype = 'acos';

%% Auto Load
LUT_depth_base = 16;
p_upper = fi([], 1, w_frac+2, w_frac);
p_comp_lower = fi([], 0, w_frac+1, w_frac);
p_atan_table = fi([], 0, M+1, M+1);
p_current_phi = fi([], 1, M+3, M+1);
ppp = {p_upper p_comp_lower p_current_phi};

load(['invc_param_k' num2str(k) '.mat'], "dd_table", "intlv_table");

castR_intlv_table = cast(intlv_table, 'like', p_upper);
dd_intlv_table = double(castR_intlv_table);

%% Construct PBT with DFS
ptree_p = [0];
pref_code = {''};
sr_pos = 1;
depths = [0, 0]; % current depth, max depth
Boundary_table = cell(1, 2^(k+1));
[ptree_p, pref_code, sr_pos, mdepth, Boundary_table, IB, IS] = search_tree(ptree_p, pref_code, sr_pos, dd_intlv_table(2:end-1), depths(1), depths(2), Boundary_table, 1, 0);
% cut blank
for ii = 1:length(Boundary_table)
    if(isempty(Boundary_table{ii}))
        break;
    end
end
Boundary_table(ii:end) = [];

%% Plot PBT
figure;
treeplot(ptree_p);

fid_pt = fopen('pt.txt', 'w');
for i = 1:length(ptree_p)-1
    fprintf(fid_pt, '%d ', ptree_p(i) );
end
    fprintf(fid_pt, '%d', ptree_p(length(ptree_p)) );
fclose(fid_pt);

% Python

%% Make Boundary Table
% Level statistics
sigma_la = zeros(1, length(Boundary_table));
tau_la = zeros(1, length(Boundary_table));
for ii = 1:length(Boundary_table)
    if(isempty(Boundary_table{ii}))
        break;
    end
    bt_temp = Boundary_table{ii};
    IS = bt_temp(1);
    sigma_la(ii) = bt_temp(2);
    tau_la(ii) = bt_temp(3);
end
tbl = tabulate(sigma_la);
offset = 0;
for ii = 1:height(tbl)
    if(tbl(ii-offset, 2) == 0)
        tbl(ii-offset, :) = [];
        offset = offset + 1;
    end
end
% Level traversal
tblc = cell(height(tbl), 1);
round_tail0_width = zeros(1, height(tbl));
round_tail0_width_ib = zeros(1, length(Boundary_table));
level_end_num = zeros(1, height(tbl));
st_num = 1;
max_head_width = 0;
for ii = 1:height(tbl)
    tblc{ii} = st_num:st_num+tbl(ii, 2)-1;

    % tail width
    mt0 = 41;
    for jj = st_num:st_num+tbl(ii, 2)-1
        bt_temp = Boundary_table{jj};
        str_bin = bin(castR_intlv_table(bt_temp(1)+1));
        nt0 = count_tail0(str_bin);
        if(nt0 < mt0)
            mt0 = nt0;
        end
    end
    round_tail0_width(ii) = mt0;
    round_tail0_width_ib(st_num:st_num+tbl(ii, 2)-1) = mt0;
    level_end_num(ii) = st_num+tbl(ii, 2)-1;

    % head width
    if(ii == height(tbl))
        break;
    end
    bt_temp = Boundary_table{st_num+tbl(ii, 2)-1};
    str_bin1 = bin(castR_intlv_table(bt_temp(1)+1));
    bt_temp = Boundary_table{st_num+tbl(ii, 2)};
    str_bin2 = bin(castR_intlv_table(bt_temp(1)+1));
    head_width = diff_head(str_bin1, str_bin2);
    if(head_width > max_head_width)
        max_head_width = head_width;
    end

    st_num = st_num+tbl(ii, 2);
end
% Leaf node traversal
tab_Is = zeros(1, length(Boundary_table));
tab_bmid = cell(1, length(Boundary_table));
str_len = zeros(1, length(Boundary_table));
fprintf('Boundary Table:\n');
fprintf('Level 1:\n');
lv_count = 1;
for ii = 1:length(Boundary_table)
    bt_temp = Boundary_table{ii};
    IS = bt_temp(1);
    sigma_l = bt_temp(2);
    tau_l = bt_temp(3);
    str_bin = bin(castR_intlv_table(IS+1));
    idx = bin2dec(str_bin(3:2+sigma_l)) - tau_l + 1; % skip sign bit
    fprintf('IB:%3d, Idx:%3d, IS:%3d, Tau:%5d  ', ii, idx, IS, tau_l);
    fprintf(['Prefix+bmid+tail0: ' str_bin(3:2+sigma_l) ' ' str_bin(3+sigma_l:end-round_tail0_width_ib(ii)) ' ' str_bin(end-round_tail0_width_ib(ii)+1:end) '  Tau(BIN)' dec2bin(tau_l) '\n']);
    if(sum(ii == level_end_num) && ii ~= length(Boundary_table))
        lv_count = lv_count + 1;
        fprintf('\n');
        fprintf('Level %d:\n', lv_count);
    end

    tab_bmid{ii} = str_bin(3+sigma_l:end-round_tail0_width_ib(ii));
    str_len(ii) = length(str_bin)-round_tail0_width_ib(ii) - (3+sigma_l) + 1;
    tab_Is(ii) = IS;
end
max_len_str = max(str_len);
fprintf('\n');
%% Write Coe
%  Note1: In Verilog, the starting index begins at 0 rather than 1.
if(~exist([arctype '_coe_files_M' num2str(M)], "dir"))
    mkdir([arctype '_coe_files_M' num2str(M)]);
end

% Boundary Table
% fid_comb_table = fopen(['./' arctype '_coe_files_M' num2str(M) '/bound_tab_k' num2str(k) '.coe'], 'w');
% fprintf(fid_comb_table, 'memory_initialization_radix=2;\nmemory_initialization_vector=');
% for i = 1:length(tab_Is)
%     str_idx = tab_Is(i);
%     if(str_idx == 0)
%         str_idx = 512;
%     end
%     if(length(tab_bmid{str_idx}) < max_len_str)
%         str_temp = [repmat('0', [1 max_len_str-length(tab_bmid{str_idx})]) tab_bmid{str_idx}];
%     else
%         str_temp = tab_bmid{str_idx};
%     end
%     fprintf(fid_comb_table, '%s,\n', [dec2bin(tab_Is(i)-1, k) str_temp]); % Note1
% end
% lut_depth_B = ceil(length(tab_Is)/LUT_depth_base)*LUT_depth_base;
% for i = length(tab_Is):lut_depth_B
%     fprintf(fid_comb_table, [repmat('0', [1 max_len_str+k]) '\n']);
% end
% fprintf(fid_comb_table, ';');
% 
% fclose(fid_comb_table);
% fprintf('LUT(Boundary Table) Config-Width x Depth: %d, %d\n', max_len_str+k, lut_depth_B);



% Leaf Node Table
% fid_comb_table = fopen(['./' arctype '_coe_files_M' num2str(M) '/Leaf_Node_Tab.coe'], 'w');
% fprintf(fid_comb_table, 'memory_initialization_radix=2;\nmemory_initialization_vector=');
% for i = 1:length(head_table)
%     str_idx = head_table(i);
%     if(str_idx == 0)
%         str_idx = 512;
%     end
%     if(length(str_table{str_idx}) < max_len_str)
%         str_temp = [repmat('0', [1 max_len_str-length(str_table{str_idx})]) str_table{str_idx}];
%     else
%         str_temp = str_table{str_idx};
%     end
%     fprintf(fid_comb_table, '%s,\n', [dec2bin(head_table(i)-1, 9) str_temp]);
% end
% fprintf(fid_comb_table, ';');


% Solution Space Table
tab_depth = length(dd_table);
tab_phi = cast(zeros(1, tab_depth), 'like', p_current_phi);
tab_upper = cast(zeros(1, tab_depth), 'like', p_upper);
tab_lower = cast(zeros(1, tab_depth), 'like', p_upper);

for ti = 1:tab_depth
    [tab_upper(ti), tab_lower(ti), tab_phi(ti)] = invc_acos_dd_actlen(dec2bin(dd_table(ti), k+1), ppp, arctype);
end
save([arctype '_invc_tab.mat'], "dd_table", "intlv_table", "tab_upper", "tab_lower", "tab_phi");

fid_SP_table = fopen(['./' arctype '_coe_files_M' num2str(M) '/Solution_Space_Tab_k' num2str(k) '.coe'], 'w');
fprintf(fid_SP_table, 'memory_initialization_radix=2;\nmemory_initialization_vector=');
for i = 1:tab_depth
    uup = bin(tab_upper(i));
    uul = bin(tab_lower(i));
    uuh = bin(tab_phi(i));
    fprintf(fid_SP_table, '%s,\n', [uup(3:end)  uul(3:end)  uuh(2:end)]);
end
lut_depth_SP = ceil(tab_depth/LUT_depth_base)*LUT_depth_base;
width_SP = p_upper.FractionLength*2+p_current_phi.WordLength-1;
for i = tab_depth+1:lut_depth_SP
    fprintf(fid_SP_table, [repmat('0', [1 width_SP]) '\n']);
end
fprintf(fid_SP_table, ';');
fclose(fid_SP_table);
fprintf('LUT(Solution Space Table) Config-Width x Depth: %d, %d\n', width_SP, lut_depth_SP);

%% Fixed Design
fprintf('Head Bit Range: [40, %d]\n', 43 - max_head_width);
head_range = zeros(height(tbl), 2);
for i = 1:height(tbl)
    list = tblc{i};
    bt_temp = Boundary_table{list(1)};
    str_bin1 = bin(castR_intlv_table(bt_temp(1)+1));
    bt_temp = Boundary_table{list(end)};
    str_bin2 = bin(castR_intlv_table(bt_temp(1)+1));
    st_head = bin2dec(str_bin1(1:max_head_width));
    ed_head = bin2dec(str_bin2(1:max_head_width));
    fprintf(['Case ' num2str(i) ' (' num2str(st_head) '-' num2str(ed_head) '):\n'])
    head_range(i, :) = [st_head ed_head];
end


%% Generate Ref Verilog
width_head = max_head_width-2;
fprintf('|-------------------- Verilog --------------------|\n');
for i = 1:height(tbl)
    temp_list = tblc{i};
    bt_temp = Boundary_table{temp_list(1)};
    str_bin1 = bin(castR_intlv_table(bt_temp(1)+1));
    bt_temp = Boundary_table{temp_list(end)};
    str_bin2 = bin(castR_intlv_table(bt_temp(1)+1));
    diff_num = diff_head(str_bin1(3:end), str_bin2(3:end));
    common_prefix = str_bin1(3:2+diff_num-1);
    width_cp = length(common_prefix);

    sigma_l = sigma_la(temp_list(1));
    tau_l = tau_la(temp_list(1));
    if(i == height(tbl))
        fprintf('        end else begin\n');
    else
        fprintf(['        end else if(x_head < ' num2str(width_head) '''d' num2str(head_range(i+1, 1)) ') begin\n']);
    end
    fprintf(['            x_prefix <= {' num2str(width_cp) '''b' common_prefix ', x_reg[' num2str(40-width_cp) ':' num2str(41-sigma_l) ']} - ' num2str(tau_l) ';\n']);
    fprintf(['            x_compare <= x_reg[' num2str(40-sigma_l) ':' num2str(round_tail0_width(i)) '];\n']);

end






function [ptree_p, pref_code, sr_pos, mdepth, Boundary_table, IB, IS] = search_tree(ptree_p, pref_code, sr_pos, intervl_array, cdepth, mdepth, Boundary_table, IB, IS)
    temp_pcode = pref_code{sr_pos};
    temp_range = prefix2range(temp_pcode);
    in_num = isin_range(temp_range, intervl_array);
    if(in_num > 1)
        parent_pos = sr_pos;
        % left child
        [ptree_p, pref_code, sr_pos, mdepth, Boundary_table, IB, IS] = search_tree([ptree_p parent_pos], [pref_code [temp_pcode '0']], sr_pos+1, intervl_array, cdepth+1, mdepth, Boundary_table, IB, IS);
        % right child
        [ptree_p, pref_code, sr_pos, mdepth, Boundary_table, IB, IS] = search_tree([ptree_p parent_pos], [pref_code [temp_pcode '1']], sr_pos+1, intervl_array, cdepth+1, mdepth, Boundary_table, IB, IS);
        mdepth = max(mdepth, cdepth);
    else
        mdepth = max(mdepth, cdepth);
        IS = IS+in_num;
        if(mdepth == cdepth)
            % Boundary Table:    Is \sigma_l \tau_l b_{mid}
            offset = bin2dec(pref_code{end})+1 - IB;
            Boundary_table{IB} = [IS cdepth offset ];
            IB = IB+1;
        else
            child_num = 2^(mdepth - cdepth);
            offset = bin2dec(pref_code{end})*child_num+1 - IB;
            for ii = IB:IB+child_num-1
                Boundary_table{ii} = [IS mdepth offset ];
            end
            IB = IB + child_num;
        end
        return;
    end
end

function [range] = prefix2range(pcode)
    range = [0 1];
    if(isempty(pcode))
        return;
    else
        ll = length(pcode);
        for i = 1:ll
            if(strcmp(pcode(i), '0'))
                range(2) = range(2) - (range(2) - range(1))/2;
            else
                range(1) = range(1) + (range(2) - range(1))/2;
            end
        end
    end
end


function [in_num] = isin_range(range, data_array)
    in_idx = (data_array >= range(1)) .* (data_array < range(2));
    in_num = sum(in_idx);
end

function [tail0_num] = count_tail0(str)
    ls = length(str);
    for i = ls:-1:1
        if(strcmp(str(i), '1'))
            tail0_num = ls - i;
            return;
        end
    end
    tail0_num = ls;
end

function [diff_num] = diff_head(str1, str2)
    ls = length(str1);
    for i = 1:ls
        if(strcmp(str1(i), str2(i)) == 0)
            diff_num = i;
            return;
        end
    end
    diff_num = 0;
end

function [upper, lower, current_phi] = invc_acos_dd_actlen(dd_idx, ppp, arctype)
    atan_table = zeros(1, 100); % 初始表
    for i = 1:100
        atan_table(i) = atan(2^(-(i-1)));
    end

    p_upper = ppp{1};
%     p_comp_lower = ppp{2};
    p_current_phi = ppp{3};

    len_idx = length(dd_idx);
    ddd = zeros(1, len_idx);
    ddd(1) = 1;
    for ii = 2:len_idx
        if(dd_idx(ii) == '1')
            ddd(ii) = 1;
        else
            ddd(ii) = -1;
        end
    end
    if(strcmp(arctype, 'acos'))
        ddd(2:end) = -ddd(2:end);
    end
    acc_phi = sum(atan_table(1:len_idx).*ddd);
    current_phi = cast(acc_phi, 'like', p_current_phi);
    lower = cast(cos(double(current_phi)), 'like', p_upper);
    upper = cast(sin(double(current_phi)), 'like', p_upper);
end