%Om Namah sivaya
%******************************************************%
%ELD using Cubic cost function with transmission losses
% for 24 time intervals and 5 runs
%*******************************************************%
% gbestmin = Power distribution (Pi).
% min_cost = Minimum cost.
% cost_run = minimum cost for each run.
% convergence = Convergence of cost of all run.
clc
clear
delete convergence1.dat interval_cost.dat load_distribution_interval.dat runno_cost.dat;
global data  Pd interval
interval = 1.0;
for i_interval = 1:interval
[popsize, MAXITER, dimension, w2, w1, runno, data, Pd, interval] = my_input();

% global data B  Pd interval


sum1 = 0;
sum2 = 0;
mean = 0;
totaltime = 0;
data1 = zeros(runno,MAXITER);

Pd(i_interval) = Pd(i_interval) ;% since value of B00 is not there in Basu et al
% Pd(i_interval) = Pd(i_interval) + B00;
l = data(:,1)';
u = data(:,2)';
lu = [l' u'];
n = length(data(:,1));
lu = lu';
range = lu(2,:) - lu(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for run = 1:runno
    mean = 0;
    T = cputime;
    x = rand(popsize,dimension,1) ;
    for k = 1:dimension
        x(:,k) = x(:,k)*range(k) + lu(1,k);
    end
    pbest = x;
    gbest = zeros(1,dimension);
    for i = 1:popsize
        %f_x(i) = f6_cubic_ELD_5unit(x(i,:),i_interval);
        [f_x(i), E_cost, F_cost ] = evalpopu(x(i,:),i_interval);
        f_pbest(i) = f_x(i);
    end
    g = min(find(f_pbest == min(f_pbest(1:popsize))));%g=find(f_pbest==min(f_pbest(1:popsize)));
    gbest = pbest(g,:);
    f_gbest = f_pbest(g);
    MINIUM = f_pbest(g);
    for t = 1:MAXITER
        disp('Minimum=');
        disp(MINIUM);
        disp('Iteration=');
        disp(t);
        disp('run no=');
        disp(run);

        out(t,run) = MINIUM; %for convergence
        out1(t,run)= F_cost;
        out2(t,run)= E_cost;
        beta = (w2-w1)*(MAXITER-t)/MAXITER+w1;
        mbest = sum(pbest)/popsize;

            cost_run(run,:) = f_gbest; %cost_run = minimum cost for each run.
            [min_cost,min_index] = min(cost_run); %min_cost = minimum cost of all runs.

             convergence1(t,:) = t;
             convergence2 = [out(:,min_index),out1(:,min_index),out2(:,min_index)];
%              convergence2 = [convergence1 out(:,min_index)];

            %convergence = [t,out(:,min_index)];%convergence of the minimum fuel cost.
            for j = 1:dimension %min_index = minimum cost index of all runs.
                data2(j,run) = gbest(j);
                gbestmin(j) = data2(j,min_index);
%                 gbestmin = gbestmin1';
            end
        
            for i = 1:popsize  

%**************************************************************************            
%**********Different probability distribution functions used************
    %        fi = rand(1,dimension);             % Uniform distribution
    %          fi = abs(randn(1,dimension));       % Normal distribution
  %        fi = abs(sum(trnd(1,dimension)));   % Cauchy distribution
     %        fi = abs(sum(exppdf(1,dimension))); % Exponential distribution
        fi = abs(sum(chi2pdf(1,dimension)));% chisquare distribution
%           fi = abs(sum(wblrnd(1,dimension))); % Weibull distribution   
%**************************************************************************           
            p = fi.*pbest(i,:) + (1-fi).*gbest;
            u = rand(1,dimension);
            b = beta*abs(mbest - x(i,:));
            v = -log(u);
            %y=p+((-1).^ceil(w1+u)).*b.*v;%
            y = p+((-1).^ceil(w1+rand(1,dimension))).*b.*v;
            x(i,:) = y;
            %x(i,:)=sign(y).*min(abs(y),xmax); 
          
             for j=1:1;
             if x(i,j)<10;
                 x(i,j)=10;
            end;
            if x(i,j)>125.00;
                x(i,j)=125.00;
            end;
             %x(i,j)=round(x(i,j));
         end;
         for j=2:2;
             if x(i,j)<10;
                 x(i,j)=10;
            end;
            if x(i,j)>150.00;
                x(i,j)=150.00;
            end;
             %x(i,j)=round(x(i,j));
         end;
         for j=3:3;
             if x(i,j)<35;
                 x(i,j)=35;
            end;
            if x(i,j)>225.00;
                x(i,j)=225.00;
            end;
             %x(i,j)=round(x(i,j));
         end;
         for j=4:4;
             if x(i,j)<35;
               
                 x(i,j)=35;
            end;
            if x(i,j)>210.00;
                x(i,j)=210.00;
            end;
             %x(i,j)=round(x(i,j));
         end;
        for j=5:5;
             if x(i,j)<130;
                x(i,j)=130;
            end;
            if x(i,j)>325.00;
                x(i,j)=325.00;
            end;
             %x(i,j)=round(x(i,j));
         end; 
         for j=6:6;
             if x(i,j)<125;
                x(i,j)=125;
            end;
            if x(i,j)>315.00;
                x(i,j)=315.00;
            end;
             %x(i,j)=round(x(i,j));
         end; 

            %f_x(i) = f6_cubic_ELD_5unit(x(i,:),i_interval);
            [f_x(i), E_cost, F_cost ] = evalpopu(x(i,:),i_interval);
            if f_x(i) < f_pbest(i)
                pbest(i,:) = x(i,:);
                f_pbest(i) = f_x(i);
            end
            if f_pbest(i) < f_gbest
                gbest = pbest(i,:);
                f_gbest = f_pbest(i);
            end            
            MINIUM = f_gbest;                    
        end
        disp('Power Output=');
        disp(gbest);
    
        data1(run,t) = MINIUM;
        if MINIUM < 2500
            mean = mean+1;
        end
    end
    sum1 = sum1 + mean;  
    sum2 = sum2 + MINIUM;
    %MINIUM
    time = cputime - T;
    totaltime = totaltime + time;

end
disp('Minimum=');
disp(MINIUM);
% dec = var(data1);   %·½²î
% St_dev = sqrt(dec);    
%gbestmin = min(data2,[],2)
Fmin = min(data1)
% St_dev;
%f_gbest
totaltime
%convergence = zeros(interval,dimension+1);
 convergence = [convergence1 convergence2]
% convergence2 = out(:,min_index);
% convergence = cat(2,[convergence1,convergence2]);
clear mean
rho = out(MAXITER,:)';
theta = (0:2*pi/29:2*pi)';
std_dev=std(rho);
means=mean(rho);
rho1 = (rho - min(rho))/(max(rho)-min(rho));
polarplot(theta,rho1,'*');
hold on
Mn=(means - min(rho))/(max(rho)-min(rho));
Mn1(1:30)=Mn;
polarplot(theta,Mn1,'--');
hold off
AA=[theta*180/pi,rho1,Mn1']
BaseName7='polar';
FileName7=[BaseName7,'.dat']
fileID7 = fopen(FileName7,'w');
fprintf(fileID7,'%16.8f\t %16.8f\t %16.8f \n', AA');



BaseName='convergence';
FileName=[BaseName,'.dat']
% fileID1 = fopen(FileName,'w');
% fprintf(fileID1,'%3d %16.8f\n', convergence');
FileName=[BaseName,num2str(i_interval),'.dat']
fileID3 = fopen(FileName,'w');
fprintf(fileID3,'%s \t%s \t\t%s \t%s \n', 'iter no','total cost', 'fuel cost', 'emission cost');
fprintf(fileID3,'%3d\t %16.8f \t%16.8f \t%16.8f \n', convergence');
% fileID1 = fopen('load.dat','w');
% fileID2 = fopen('total_cost.dat','w');
% fprintf(fileID1,'%2d %16.8f %16.8f %16.8f %16.8f\n', B');
fileID2 = fopen('load_distribution_interval.dat','a+');
fprintf(fileID2,'%2d \t %16.8f\t %16.8f\t %16.8f\t %16.8f\t %16.8f\t %16.8f\t  %16.8f\n', i_interval,gbestmin',MINIUM);
fileID5 = fopen('interval_cost.dat','a+');
fprintf(fileID5,'%s \t%s \t\t\t%s \t\t%s\n', 'interval','interval_cost','std_dev','means');
fprintf(fileID5,'%2d \t\t%16.8f \t%16.8f \t%16.8f\n', i_interval,convergence2(MAXITER),std_dev,means);
% dlmwrite('convergence_.dat',cat(2,[convergence1,convergence2]),'Delimiter','\t')
end
fileID4 = fopen('runno_cost.dat','a+');
fprintf(fileID4,'%16.8f\n', cost_run');