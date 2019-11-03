%% 该函数演示多目标perota优化问题
%清空环境
clc
clear
dataset=xlsread('dataset.xlsx');
load data
%% 初始参数
%objnum=size(P,1); %类中物品个数
%weight=92;        %总重量限制

%初始化程序
Dim=10;     %粒子维数    路径有10个服务组成
xSize=50;  %种群个数       
MaxIt=200; %迭代次数
c1=0.8;    %算法参数
c2=0.8;    %算法参数 
wmax=1.2;  %惯性因子
wmin=0.1;  %惯性因子

x=unidrnd(1000,xSize,Dim);  %粒子初始化
v=zeros(xSize,Dim);      %速度初始化

xbest=x;           %个体最佳值
gbest=x(1,:);      %粒子群最佳位置

% 粒子适应度值 
cost=zeros(1,xSize);   %粒子成本目标
time=zeros(1,xSize);   %粒子时间目标
reliable=zeros(1,xSize);   %可靠性
ec=zeros(1,xSize);   %能耗
% 最优值初始化
cobest=zeros(1,xSize); %粒子成本目标
tibest=zeros(1,xSize); %粒子时间目标
rebest=zeros(1,xSize);  %可靠性
ecbest=zeros(1,xSize);  %能耗 ++++++++


% 上一次的值
coPrior=zeros(1,xSize);%粒子价值目标
tiPrior=zeros(1,xSize);%粒子体积目标
rePrior=zeros(1,xSize);%记录重量，以求约束
ecPrior=zeros(1,xSize);  %记录重量，以求约束+++++++

%计算初始目标向量
for i=1:xSize
    for j=1:Dim %控制类别
        cost(i) = cost(i)+dataset(x(i,j),1);  %粒子成本
        time(i) = time(i)+dataset(x(i,j),2);  %粒子时间
        reliable(i) = reliable(i)+dataset(x(i,j),3);  %粒子可靠性
        ec(i) = ec(i)+dataset(x(i,j),4);  %粒子能耗+++++++
    end
end
% 粒子最优位置
cobest=cost;tibest=time;rebest=reliable;ecbest=ec;

%% 初始筛选非劣解
flj=[];
fljx=[];
fljNum=0;
%两个实数相等精度
tol=1e-7;
for i=1:xSize
    flag=0;  %支配标志
    for j=1:xSize  
        if j~=i
            if ((cost(i)>cost(j)) &&  (time(i)>time(j)) && reliable(i)<reliable(i) && ec(i)>ec(i)) 
                flag=1;%不支配标志
                break;
            end
        end
    end
    
    %判断有无被支配
    if flag==0
        fljNum=fljNum+1;
        % 记录非劣解
        flj(fljNum,1)=cost(i);flj(fljNum,2)=time(i);flj(fljNum,3)=reliable(i);flj(fljNum,4)=ec(i);
        % 非劣解位置
        fljx(fljNum,:)=x(i,:); 
    end
end

%% 循环迭代
for iter=1:MaxIt
    
    % 权值更新
    w=wmax-(wmax-wmin)*iter/MaxIt;
     
    %从非劣解中选择粒子作为全局最优解
    s=size(fljx,1);       
    index=randi(s,1,1);  
    gbest=fljx(index,:);

    %% 群体更新
    for i=1:xSize
        %速度更新
        v(i,:)=w*v(i,:)+c1*rand(1,1)*(xbest(i,:)-x(i,:))+c2*rand(1,1)*(gbest-x(i,:));
        
        %位置更新
        x(i,:)=x(i,:)+v(i,:);  
        index1=find(x(i,:)<=0);
        if ~isempty(index1)
            x(i,index1)=rand(size(index1));
        end
        for j=1:10
            if x(i,j)>1000
                x(i,j)=1000;
            end
            
        x(i,:)=ceil(x(i,:));        
        end
    end
    
    %% 计算个体适应度
    coPrior(:)=0;
    tiPrior(:)=0;
    rePrior(:)=0;
    ecPrior(:)=0;
    for i=1:xSize
        for j=1:Dim %控制类别
            coPrior(i) = coPrior(i)+dataset(x(i,j),1);  %计算粒子成本
            tiPrior(i) = tiPrior(i)+dataset(x(i,j),2);  %计算粒子时间
            rePrior(i) = rePrior(i)+dataset(x(i,j),3);  %计算粒子可靠性
            ecPrior(i) = ecPrior(i)+dataset(x(i,j),4);  %计算粒子能耗
        end
    end
    
    %% 更新粒子历史最佳
    for i=1:xSize
        for j=1:xSize
        %现在的支配原有的，替代原有的
         if ((cost(i)<cost(j)) &&  (time(i)<time(j)) && reliable(i)>reliable(j) && ec(i)<ec(j))
                xbest(i,:)=x(i,:);%没有记录目标值
                cobest(i)=coPrior(i);tibest(i)=tiPrior(i);rebest(i)=rePrior(i);ecbest(i)=ecPrior(i);
        elseif rand(1,1)<0.5
                xbest(i,:)=x(i,:);
                  cobest(i)=coPrior(i);tibest(i)=tiPrior(i);rebest(i)=rePrior(i);ecbest(i)=ecPrior(i);
         end
        end
    end

    %% 更新非劣解集合
    cost=coPrior;
    time=tiPrior;
    reliable=rePrior;
    ecliable=ecPrior;
    %更新升级非劣解集合
    s=size(flj,1);%目前非劣解集合中元素个数
   
    %先将非劣解集合和xbest合并
    cccx=zeros(1,s+xSize);
    tttx=zeros(1,s+xSize);
    rrrx=zeros(1,s+xSize);
    eeex=zeros(1,s+xSize);
    
    cccx(1:xSize)=cobest;cccx(xSize+1:end)=flj(:,1)';
    tttx(1:xSize)=tibest;tttx(xSize+1:end)=flj(:,2)';
    rrrx(1:xSize)=rebest;rrrx(xSize+1:end)=flj(:,3)';
    eeex(1:xSize)=ecbest;eeex(xSize+1:end)=flj(:,4)';
    
    xxbest=zeros(s+xSize,Dim);
    xxbest(1:xSize,:)=xbest;
    xxbest(xSize+1:end,:)=fljx;
   
    %筛选非劣解
    flj=[];
    fljx=[];
    k=0;
    tol=1e-7;
    for i=1:xSize+s
        flag=0;%没有被支配
        %判断该点是否非劣
        for j=1:xSize+s 
            if j~=i
                if ((cccx(i)>cccx(j)) &&  (tttx(i)>tttx(j)) &&  (rrrx(i)<rrrx(j)) && (eeex(i)<eeex(j)))  %有一次被支配
                    flag=1;
                    break;
                end
            end
        end

        %判断有无被支配
        if flag==0
            k=k+1;
            flj(k,1)=cccx(i);flj(k,2)=tttx(i);flj(k,3)=rrrx(i);flj(k,4)=eeex(i);%记录非劣解
            fljx(k,:)=xxbest(i,:);%非劣解位置
        end
    end
    
    %去掉重复粒子
    repflag=0;   %重复标志
    k=1;         %不同非劣解粒子数
    flj2=[];     %存储不同非劣解
    fljx2=[];    %存储不同非劣解粒子位置
    flj2(k,:)=flj(1,:);
    fljx2(k,:)=fljx(1,:);
    for j=2:size(flj,1)
        repflag=0;  %重复标志
        for i=1:size(flj2,1)
            result=(fljx(j,:)==fljx2(i,:));
            if length(find(result==1))==Dim
                repflag=1;%有重复
            end
        end
        %粒子不同，存储
        if repflag==0 
            k=k+1;
            flj2(k,:)=flj(j,:);
            fljx2(k,:)=fljx(j,:);
        end
        
    end
    
    %非劣解更新
    flj=flj2;
    fljx=fljx2;
 disp(['In iteration ' num2str(iter) ]);
    save results
end