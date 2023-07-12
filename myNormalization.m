function y = myNormalization(x, type)
% ==========================================
% function y = myNormalization(x, type)
% type = 0, 除以总的最大值
% type = 1, 各个波段除以其最大值
% ==========================================
x = double(x);
[m,n,p] = size(x);
y = zeros(m,n,p);
if type == 0
    y = x/max(x(:));
elseif type == 1
    for i = 1:p
        temp = x(:,:,i);
        y(:,:,i) = temp/max(temp(:));
    end
elseif type == 2
    for i=1:p
        otmp = x(:,:,i); 
        min_x(i) = min(otmp(:));
        max_x(i) = max(otmp(:));
        otmp = otmp - min_x(i);
        scale(i) = max_x(i)-min_x(i);
        %scale to [0,1]
        y(:,:,i) = otmp/scale(i);
    end
end