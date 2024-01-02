clear all

% time data
fid = fopen('ackermannTimeData.txt');
line = fgetl(fid);
t = [];
while line ~= -1
    line_val = sscanf(line, '%f')';
    t = [t, line_val];
    line = fgetl(fid);
end
fclose(fid);

% state data
fid = fopen('ackermannStateData.txt');
line = fgetl(fid);
x = [];
while line ~= -1
    line_val = sscanf(line, '%f')';
    x = [x; line_val];
    line = fgetl(fid);
end
x = x';
fclose(fid);

% input data
fid = fopen('ackermannInputData.txt');
line = fgetl(fid);
u = [];
while line ~= -1
    line_val = sscanf(line, '%f')';
    u = [u; line_val];
    line = fgetl(fid);
end
u = u';
fclose(fid);

% reference
fid = fopen('ackermannRefData.txt');
line = fgetl(fid);
x_ref = [];
cell_cnt = 1;
while isempty(line) || ~all(line == -1)
    line_val = sscanf(line, '%f')';
    if ~isempty(line_val)
        x_ref = [x_ref; line_val];
    else
        x_ref_cell{cell_cnt} = x_ref';
        x_ref = [];
        cell_cnt = cell_cnt+1;
    end
    line = fgetl(fid);
end
fclose(fid);

% first reference
x_ref0 = zeros(size(x)); % init
for i = 1:length(x)
    x_ref0(:,i) = x_ref_cell{i}(:,1);
end

% put into struct
ackermannData.t = t;
ackermannData.x = x;
ackermannData.u = u;
ackermannData.x_ref = x_ref_cell;
ackermannData.x_ref0 = x_ref0;

save ackermannData ackermannData