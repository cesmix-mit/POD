n = 1000;
m = 1000;

if m>=n
    b = 1;
    c = [0 n];
else
    b = ceil(n/m); 
    c = round(linspace(0, 1000, b+1));    
end


% a = 1;
% b = 1000;
% m = 150;
% n = (b-a+1)/m;
%m = (b-a)/(n-1);
%double delta = (end - start) / (num - 1);

