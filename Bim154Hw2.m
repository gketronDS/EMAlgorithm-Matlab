%% 5c
r = binornd(10,0.5,10,1)
ravg = mean(r)
figure;
hist1=histogram(r);
%% 5d
expset = datasample(r,10);
expsetavg = mean(expset);
figure;
hist2 = histogram(expset);
%% 5e
x = [];
for index = 1:1000
    omniset = datasample(r,10);
    omnisetavg = mean(omniset);
    x(end+1) = omnisetavg;
end
x;
xavg = mean(x)
figure 
hist3 = histogram(x);
% ravg is very close to xavg (1000 times). The difference is only 0.0116.
%% 5f
y = [];
for index = 1:100000
    omniset = datasample(r,10);
    omnisetavg = mean(omniset);
    y(end+1) = omnisetavg;
end
yavg = mean(y)
figure 
hist4 = histogram(y);
% part f is much closer to the actual, with a difference of only 0.0025.
% part d is farther from the actual with a difference of 0.0225. The more
% samples, the more accurate. 
%% 5g
bigr = binornd(10,0.5,10000,1);
bigravg = mean(bigr)
figure 
hist5=histogram(bigr);
%% 5h
bigset = datasample(bigr,10000);
bigsetavg = mean(bigset)
figure 
hist6 = histogram(bigset);
%% 5i
z = [];
for index = 1:1000
    megaset = datasample(bigr,10000);
    megasetavg = mean(megaset);
    z(end+1) = megasetavg;
end
z;
zavg = mean(z)
figure 
hist7 = histogram(z);
% By increasing the dataset size, the samples were able to get more
% accurate estimations of the original data set. 