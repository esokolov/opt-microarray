I = 1;
a = 2;
c = 5;
f = @(b) log(I) -1 + (a * c) / ( I  * (1 + c * b)) - log((a * c) / (1 + c * b));
x = 0:0.01:10;
y = x;
for i = 1:length(x)
    y(i) = f(x(i));
end
plot(x, y)
hold on
fl = @(b, b0) log(I) - 1 + (a*c) / (I * (1 + c * b)) - log((a*c) / (1 + c * b0)) + c * (b - b0) / (1 + c * b0);
%%
yl = x;
for i = 1:length(x)
yl(i) = fl(x(i), 1.9);
end
plot(x, yl, 'r')

%%
f1 = @(b, b0) 0.5 * b .^2 * c / (1 + b0 * c) - b0 * b * a * c .^ 2 / (I * ((1 + b0 * c) .^ 2) ) + 5;
y1 = x;
for i = 1:length(x)
    y1(i) = f1(x(i), 2.5);
end
plot(x, y1, 'r')


%%
fc = @(b) (a * c) / (1 + c * b);
yc = x;
for i = 1:length(x)
    yc(i) = fc(x(i));
end
plot(x, yc, 'r')
%%
fcc = @(b) -1 - log((a * c) / (1 + c * b));
ycc = x;
for i = 1:length(x)
    ycc(i) = fcc(x(i));
end
plot(x, ycc, 'g')