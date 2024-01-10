function [fcn, E_cost, F_cost ] = evalpopu(pop1,i_interval)
%===============================================================
%  Objective function is estimated here
%===============================================================
% Coded By Dr V. Ravikumar Pandi and Ms Angel T S
% Amrita Vishwa Vidyapeetham, Kollam
% FDP on Power Quality and Smart Grid Optimization  
% Held at Bishop Jerome Institute, Kollam. 
% Date : 19-07-2018
%===============================================================

%data1=[ Pmin   Pmax    a           b        c        d        e     alpha     beta      gamma       eetta      delta]; 
 
data1=[   10 	125	    0.1525	    38.540		756.800      13.860     0.3300     0.0042;
         10	    150	    0.1060	    46.160		451.325      13.860     0.3300     0.0042;    
         35	    225	    0.0280	    40.400		1050.000     40.267    -0.5455     0.0068;   
         35 	210	    0.0355	    38.310		1243.530     40.267    -0.5455     0.0068;     
         130	325	    0.0211	    36.328		1658.570     42.900    -0.5112     0.0046;     
         125	315	    0.0180	    38.270		1356.660     42.900    -0.5112     0.0046;];
    
%        180    350     0.004531     7.3968   643.24;     %thermal power plant1  
%        180    350     0.004683     7.5629   666.27;     %thermal power plant2
%        180    350     0.004683     7.5629   666.27;     %thermal power plant3
%        180    350     0.004708     7.4767   672.77;];   %thermal power plant4 
            
    
% % Load Demand
PD=1200;
%PD=[487.50 460.00 461.00 450.00 454.25 470.00 488.00 613.75 668.00 669.75 654.50 ...
    %690.75 708.75 708.75 627.00 659.50 804.25 875.00 831.25 823.25 786.50 717.00 603.75 547.50]';
% PD=[1036 1110 1258 1406 1480 1628 1702 1776 1924 2022 2106 2150 ...
%     2072 1924 1776 1554 1480 1628 1776 1972 1924 1628 1332 1184 ]';
% Transmission loss coefficient
% coeff=[0.000049 0.000014  0.000015  0.000015 0.000016 0.000017 0.000017 ...
%     0.000018  0.000019  0.000020; 0.000014 0.000045 0.000016 0.000016 ...
%     0.000017  0.000015  0.000015 0.000016 0.000018 0.000018; 0.000015 ...
%     0.000016  0.000039  0.000010 0.000012 0.000012 0.000014 0.000014 ...
%     0.000016  0.000016; 0.000015 0.000016 0.000010 0.000040 0.000014 ...
%     0.000010  0.000011  0.000012 0.000014 0.000015; 0.000016 0.000017 ...
%     0.000012  0.000014  0.000035 0.000011 0.000013 0.000013 0.000015 ...
%     0.000016; 0.000017  0.000015 0.000012 0.000010 0.000011 0.000036 ...
%     0.000012 0.000012 0.000014 0.000015; 0.000017 0.000015 0.000014 ...
%     0.000011 0.000013 0.000012 0.000038  0.000016 0.000016 0.000018; ...
%     0.000018 0.000016 0.000014 0.000012  0.000013 0.000012 0.000016 ...
%     0.000040 0.000015 0.000016; 0.000019 0.000018 0.000016 0.000014 ...
%     0.000015 0.000014 0.000016 0.000015  0.000042 0.000019; 0.000020 ...
%     0.000018 0.000016 0.000015 0.000016  0.000015 0.000018 0.000016 ...
%     0.000019 0.000044;];

% To find h (max penalty factor)
F_cost_max = sum(data1(:,3).*data1(:,2).^2+data1(:,4).*data1(:,2)+data1(:,5));
E_cost_max = sum(data1(:,6)+data1(:,7).*data1(:,2)+data1(:,8).*data1(:,2).^2); 
        h =  F_cost_max /E_cost_max;
%PD=668;
%in=1;
%in;
% A=1000;
% B=1;

[pop_n  var_n]= size(pop1);
fitness = zeros(pop_n, 1);
for count = 1:pop_n
    x1= pop1(count,:)';
    %objective Function
% worst practice
%     cost=data1(1,3)*x1(1)*x1(1)+data1(1,4)*x1(1)+data1(1,5)+...
%          data1(2,3)*x1(2)*x1(2)+data1(2,4)*x1(2)+data1(2,5)+...
%          data1(3,3)*x1(3)*x1(3)+data1(3,4)*x1(3)+data1(3,5);
% best practice (generalization)
    F_cost = sum(data1(:,3).*x1.^2+data1(:,4).*x1+data1(:,5));%fuel cost fn with valve point effect
   E_cost = sum(data1(:,6)+data1(:,7).*x1+data1(:,8).*x1.^2);
 %=========================================================================
 %To evaluate the emission control cost factor "g"
 %============================================================================
%  for i=1:var_n
% F_cost_max(i) = sum(data1(i,3).*data1(i,2).^2+data1(i,4).*data1(i,2)+data1(i,5)+...
%             abs(data1(i,6).*sin(data1(i,7).*(data1(i,1)- (data1(i,2))))));
% E_cost_max(i)= sum(data1(i,8)+data1(i,9).*data1(i,2)+data1(i,10).*data1(i,2).^2+...
%             data1(i,11).*exp(data1(i,12).*(data1(i,2)))); 
%         h(i) =  F_cost_max(i) /E_cost_max(i);
%  end
 %==============================================================================
    % To find emission cost
 cost_emission = h * E_cost;
 %=================================================================================
    % To convert to the single objective function
     w1=0.5; %weighing factor for each obj fn
   cost = w1*F_cost' +(1-w1)*cost_emission';
 %=========================================================================== 
%  % To find transmission loss Ploss
%  for l=1:var_n
%      for r=1:var_n
%          P_loss = sum(x1(l,1)*coeff(l,r)*x1(r,1));
%      end
%  end
 %==========================================================================
    %start_up = sum (data1 (:,6)*exp())
%     penalty function is evaluated here
%    Ppb=(sum(x1)-PD(in))^2;
  %  Ppb=abs(sum(x1')-PD(i_interval)- P_loss')';
 Ppb=abs(sum(x1)-PD(i_interval));
    fcn=cost+100*Ppb;
    %fitness(count)=A/(B+fcn);
    stored_var(count,:)=[fcn  fitness(count)];
    E_cost(count,:)= [ E_cost];
    F_cost(count,:)= [ F_cost];
end
%in=in+1;
