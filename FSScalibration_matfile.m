% Script to automatically process FSS testing data to calibrate the sensor.
% This version is for test data combined into 1 mat or csv file 

% Author: Sarthak Srivastava

%%------ Input file ------------
% assumes data is loaded in the workspace as data_a matrix

% reading data from data_a 3D mat structure
[row,col,heig] = size(data_a);
data = reshape(data_a,[row*col,heig]);  %reshape to 2d
% remove all zero rows and columns
data( ~any(data,2), : ) = []; 
data( :, ~any(data,1) ) = [];

% separate into FSS1 and FSS2 and in format [I1 I2 I3 I4 beta alpha]
% (ignoring temperature and FSS_ref):
FSS1data = [data(:,3:6) data(:,1:2)];
FSS2data = [data(:,8:11) data(:,1:2)];

%FSS1error = ProcessFSSreadings(FSS1data);
FSS2error = ProcessFSSreadings(FSS2data);

% OR if csv file was provided:
%data = readmatrix('FSSDataCleanedRadians.csv'); % format [I1 I2 I3 I4 beta alpha]

%% Averaging and extracting data points
%data(:,5:6) = rad2deg(data(:,5:6)); % converting angles from rad to deg

% % There are 10 readings with the same alpha and beta values -> take
% % average and reduce the dataset
% [rows, cols] = size(data);
% i = 1; j = 1;
% while i <= rows
% fssdata(j,:) = mean(data(i:i+9,:),1); % not optimal, but convenient if no. of alpha values is not readily known (remembered!)
% j = j+1;
% i = i+10;
% end

function[ang_error] = ProcessFSSreadings(fssdata)

% INPUT: sensor data as [I1 I2 I3 I4 beta_deg alpha_deg]
% data must be unique, i.e. no duplication of points

fssdata = sortrows(fssdata,6);    % sort in ascending order of alphas

% fssdata only contains unique alpha-beta pairs -> duplicates are averaged out    
I1 = fssdata(:,1); % current values for each data point 
I2 = fssdata(:,2);
I3 = fssdata(:,3);
I4 = fssdata(:,4);
beta = fssdata(:,5); % beta angles for each data point
alpha = fssdata(:,6); % alpha angles for each data point   
no_alphas = length(unique(fssdata(:,6))); % number of alpha values in the sweep

%% Computing calibration coefficients --------------------

% Get x and y (times 2/L and 2/L respectively, ignored since constants do not matter)
% based on equations in Hamamatsu PSD datasheet
 x = ((I2 + I3) - (I1 + I4))./(I1+I2+I3+I4);
 y = ((I2 + I4) - (I1 + I3))./(I1+I2+I3+I4);   

% Compute array terms to be used later (each of these is n x 1)
x2 = x.^2;
x3 = x.^3;
x4 = x.^4;
x5 = x.^5;
x6 = x.^6;
y2 = y.^2;
y3 = y.^3;
y4 = y.^4;
y5 = y.^5;
y6 = y.^6;
xy = x.*y;
x2y = x2.*y;
xy2 = x.*y2;
xy3 = x.*y3;
x3y = x3.*y;
x2y2 = x2.*y2;
x2y3 = x2.*y3;
x3y2 = x3.*y2;
x3y3 = x3.*y3;
x4y = x4.*y;
xy4 = x.*y4;
x5y = x5.*y;
xy5 = x.*y5;
x4y2 = x4.*y2;
x2y4 = x2.*y4;

xx = -tand(alpha);  %ideal xx and yy (true values based on turntable alpha and beta)
yy = -tand(beta);

n = length(x);
      
%% Using Batch Least Squares equation (5th order)
 H = [ones(n,1) x y x2 y2 xy x3 y3 x2y xy2 x4 xy3 x3y y4 x5 x3y2 x2y3 y5 xy4 x4y];
 a = inv(H'*H)*H'*xx;  
 b = inv(H'*H)*H'*yy;
 
 FSS_xx = a(1) + a(2).*x + a(3).*y + a(4).*x2 + a(5).*y2 + a(6).*xy + a(7).*x3 + a(8).*y3 + a(9).*x2y + a(10).*xy2 + a(11).*x4 + a(12).*xy3 ...
         + a(13).*x3y + a(14).*y4 + a(15).*x5 + a(16).*x3y2 + a(17).*x2y3 + a(18).*y5 + a(19).*xy4 + a(20).*x4y ;
 FSS_yy = b(1) + b(2).*x + b(3).*y + b(4).*x2 + b(5).*y2 + b(6).*xy + b(7).*x3 + b(8).*y3 + b(9).*x2y + b(10).*xy2 + b(11).*x4 + b(12).*xy3 ...
         + b(13).*x3y + b(14).*y4 + b(15).*x5 + b(16).*x3y2 + b(17).*x2y3 + b(18).*y5 + b(19).*xy4 + b(20).*x4y  ;

    
  %% Sun vector calculation ------------
  sun_b_true = zeros(n,3);
  sun_b_measured = zeros(n,3);
  ang_error = zeros(n,1);
  
  for i = 1:n
    sun_b_true(i,:) = [xx(i), yy(i), 1]./norm([xx(i), yy(i), 1]);
    sun_b_measured(i,:) = [FSS_xx(i), FSS_yy(i), 1]./norm([FSS_xx(i), FSS_yy(i), 1]);
  
    % Sun vector measurement error:
    ang_error(i) = acosd(dot(sun_b_true(i,:),sun_b_measured(i,:))); %in degrees 
  end
  
    %% Plotting -------------------------
   [Alph,Bet] = meshgrid(-60:1:60);
   Err =griddata(alpha,beta,ang_error,Alph,Bet);
   figure
   surf(Alph,Bet,Err)
   xlabel('alpha')
   ylabel('beta')
   zlabel('angular error (deg)')
   title('FSS angular error in degrees')
   
   figure
   contourf(Alph,Bet,Err,50,'EdgeColor','None')
   colorbar
   xlabel('alpha (deg)')
   ylabel('beta (deg)')
   title('FSS angular error in degrees')
end
  