%2016_07_12增加偏移量的func_ga_06
addpath('E:\毕业\代码\函数\');
clear;
close all;
dirname='E:\毕业\数据\数据2017-04-17\';
filename=dir(dirname);
% z=load('E:\毕业\数据\gaoyongli-2016-10-18-02.xls');
for p=1:length(filename)
    if exist(filename(p).name,'file')
        continue;
    end
    z=load([dirname,filename(p).name]);
    [a,b]=size(z);
    if b==2
        az=69.25*z(:,3)-19.74;
    else
        az=69.25*z(:,3)-19.74;
    end
    for i=1:length(az)
        if isnan(az(i))==1
            az(i)=0;
        end
    end
    az(az==0)=[];

    params=[];
    [sm_03,sm_05,sm_06,segerrpo_01,segerrpo_02,Onset_loc,SigPeak_loc]=func_getadapterline(az);

    figure;
    plot(sm_03,sm_05);hold on;plot(sm_03(SigPeak_loc),sm_05(SigPeak_loc),'o');
    set(gca,'XDir','reverse')%对X方向反转;

    figure;
    plot(sm_03(Onset_loc),sm_05(SigPeak_loc)-sm_05(Onset_loc),'o--');
    set(gca,'XDir','reverse')%对X方向反转;

    figure;
    plot(sm_03(SigPeak_loc),sm_05(SigPeak_loc),'o--');
    set(gca,'XDir','reverse')%对X方向反转;
    figure;
    plot(sm_03(Onset_loc),sm_05(Onset_loc),'o--');
    set(gca,'XDir','reverse')%对X方向反转;


    %%提取包络函数：包括一些奇异点去除代码：
    [Env,Env1,Env2]=func_wave2env_02(sm_03,sm_05,segerrpo_01,segerrpo_02,Onset_loc,SigPeak_loc);

    figure;
    plot(Env(:,1),Env(:,2),'o--');set(gca,'XDir','reverse')%对X方向反转;
    figure;
    plot(Env1(:,1),Env1(:,2),'o--');set(gca,'XDir','reverse')%对X方向反转;
    figure;
    plot(Env2(:,1),Env2(:,2),'o--');set(gca,'XDir','reverse')%对X方向反转;



    % 有些数据具有较强的呼吸干扰，这个时候需要其他一些方法去除（如果用基于特征点的方法）

    %插值:140-40，50个点
    [env_interp]=func_interp(Env);
    [env_interp1]=func_interp(Env1);
    [env_interp2]=func_interp(Env2);

     env_interp(isnan(env_interp(:,2)),:)=[];




    %butter低通滤波，去除干扰：
    [n,Wn] = buttord(0.2/1,0.8/1,3,60);
    [Bb,Ba] = butter(n,Wn);
    % env_interp(:,2)=awgn(env_interp(:,2),10*log10(500));%10%
    Bf=filtfilt(Bb,Ba,env_interp(:,2)); % 进行低通滤波
    Bf=Bf;


    global Env_last 


      Env_last=[env_interp(:,1) Bf/max(Bf)];
    %  Env_last=[env_interp(:,1) env_interp(:,2)/max(env_interp(:,2))];
    %这里我们没有对数据进行平滑：


    %% 遗传算法的计算
    len=1;
    disp('遗传算法');
    resu=zeros(6,len);g=zeros(1,len);c=zeros(1,len);
    g=zeros(1,len);
    Env=Env_last;
    for v=1:5
    %     Env_last(:,2)=awgn(Env_last(:,2),10*log10(500));%5%
        Env_last=Env(int16(1+(v-1)*int16((find(Env(:,2)==max(Env(:,2)))-1)/5)):end,:);
        for i=1:len
            [resu(:,i),g(i),c(i),d{i},e{i},f{i}]=func_ga_11(6,[ 80 50 0.01,0.001 0.02 0.2],[ 200 130 0.3 0.2 0.5 1],200,1600,1e-8,1e-5);
            resu(:,i);
            g(i);
        end

        loc_seed=find(g==min(g));
        loc_seed=loc_seed(1);
        resu(:,loc_seed)
        min(g)
        SBP=resu(1,loc_seed);
        DBP=resu(2,loc_seed);
        a=resu(3,loc_seed);
        b=resu(4,loc_seed);
        A=resu(5,loc_seed);
        A2=resu(6,loc_seed);
        params=[params resu(:,loc_seed)];
        %%
        y_=[];%% 上次数据清除
        for j=1:length(Env_last(:,1))
            if SBP<Env_last(j,1)
                y_(j)=A2*(exp(a*(SBP-Env_last(j,1)))-exp(a*(DBP-Env_last(j,1)))+A)*(Env_last(j,1)+760)/300;
            elseif SBP>Env_last(j,1) && DBP<Env_last(j,1)
                y_(j)=A2*(1+a/b*(1-exp(-b*(SBP-Env_last(j,1))))-exp(a*(DBP-Env_last(j,1)))+A)*(Env_last(j,1)+760)/300;
            elseif DBP>Env_last(j,1)
                y_(j)=A2*(a/b*(1-exp(-b*(SBP-Env_last(j,1))))-a/b*(1-exp(-b*(DBP-Env_last(j,1))))+A)*(Env_last(j,1)+760)/300;
            end
        end

        figure;plot(Env_last(:,1),max(Bf)*y_);

        hold on;plot(Bf(:,1),Bf,'o--');

        hold on;plot(env_interp(:,1),env_interp(:,2),'o--');
        set(gca,'XDir','reverse')%对X方向反转;
        close all;
 %%     
        dirpic='../picture/';
        mkdir(dirpic);
        for num=1:8
             saveas(figure(num),[dirpic,filename(p).name,'+',num2str(num),'.fig']);
             saveas(figure(num),[dirpic,filename(p).name,'+',num2str(num),'.tiff']);
        end
        
    end
end
