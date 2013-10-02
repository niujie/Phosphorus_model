function p = rainfall(t)

p = 0.0;
date = datestr(t + datenum('01-Jan-1998'));
month = date(4:6);
switch month
    case {'May', 'Jun', 'Jul', 'Aug', 'Sep'}
        times = poissrnd(0.46);
        for i = 1 : times
            p = p + exprnd(1.26);
        end
%         pd1 = makedist('Poisson', 'lambda', 0.46);
%         pd2 = makedist('Exponential', 'mu', 1.26);
    otherwise
        times = poissrnd(0.05);
        for i = 1 : times
            p = p + exprnd(0.4);
        end
%         pd1 = makedist('Poisson', 'lambda', 0.05);
%         pd2 = makedist('Exponential', 'mu', 0.4);
end
% times = random(pd1);
% for i = 1 : times
%     p = p + random(pd2);
% end