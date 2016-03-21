function [att,risk,gmpa,riskd]=sha(inputfile, trl,sigma,byes)
%  SHA  seismic hazard analysis code for CPSHA
% 2015/08/25： 修改了内部算法，不再保存中间结果到文本文件.txt
% 显著提高了计算速度；是否需要再进一步？先计算后校正
%{
 INPUT:
 inputfile         ----  input file, include seismic province/source/rate/spational distribution function/Mu...
 trl                 ----  truncate level
 sigma           ----  standard devation of ln(Sa)
 byes             ----  whether consider background contribution, 1: yes;  0, no

OUTPUT:
 sha_dis.mat   ----  
 risk       ---   衰减关系不确定性校正后的地震动超越概率，对应于att
 riskd          ---    给定的超越概率
 gmpa           ---    地震动参数的幅值 
 risk_median    ---  中值地震动衰减关系计算出的超越概率
 att            ---   预定的地震动幅值
 
 DEFINE VARIABLES: 
 trl                   截断水平    
 sigma            地震动衰减关系标准差,自然对数
 deps(30)         地震带震源平均深度，30为地震带个数
 amu(30)          地震带的震级上限
 beta(30)         地震带的b值
 nz               地震带个数 number of  zones
 nz1(30)           地震带nz内潜源编号的最大值 
 rate(30)          地震带的地震年平均发生率v4(整个地震带的统计结果，不除以面积)
 hks(2,2100)     长短轴方向矩阵，2100为潜源个数,(衰减关系长短轴方向，不同潜源可以不同，
                        以中国的实践来看，长轴通常为断裂带的走向或宏观震中密集带)
 phs(2,2100)     长短轴方向概率矩阵 
 ae(2100)         潜源面积(km^2)
 ord                 对当前场点is有影响的潜源编号
 ord2(2100)      衰减关系指标，即潜源对应的衰减关系（1东部ca,以经度平均是否大于105判断) 
 xs(12,2100)     潜源顶点横坐标（经度），12为潜源顶点个数最大值
 ys(12,2100)     潜源顶点纵坐标（纬度 (潜源都为以宏观震中控制的面源，最多为12边形)
 ammas(2100)     潜源震级上限(ms)
 q0(8,2100)      潜源空间分布函数，8为震级档个数
 amstep          实际计算震级档大小 0.2ms
 ek              计算地震发生频率时的分母，即根据g-r关系，积分时使概率为1的常数
 brs             背景源对场点的影响
 ra0/rb0         衰减关系中的r0=ca8*exp(ca9*m)
 s(j,jj)         背景源影响场面积，当前震级档震级jj，当前地震动j
 brtate          背景源单位面积地震年平均发生率
 lna=ca1+ca21*m+ca22*m*m-ca3*ln(r+ca8*exp(ca9*m))+epsilon*std   ca长轴；cb短轴
 siga(4)         衰减关系对数标准差，自然对数
 ff(i,j)         潜源i中,震级档(实际计算)j单位面积上的地震年平均发生率
 amq(30)         每个震级档的起始震级,例如：4.0 5.5 6.0 6.5 7.0 7.5
 ord(2100)       对场点地震动有影响的潜源的编号，潜源统一编号，每个场点的ord可变，潜源编号从小到大排列，例如：3, 28, 78, 100, 35……
 sitex，sitey    场点经纬度坐标 ，
 b               潜源顶点pi, pi+1组成直线在y轴，垂直破裂方向的截距， 
 q0              空间分布函数，最多10个震级档
 risk            场点地震危险性矩阵，40个地震动幅值，2100个潜源
 g(20)           潜源顶点坐标wx,wy所组成直线的正切值
 rs(45)          衰减关系不确定性校正的中间变量
 ax(40)          衰减关系不确定性校正的中间变量
 brs(40)         background source贡献
 lst             潜源类型，1：断层源，2：二类面源，3，三类面源
 x，y            当前潜源顶点的坐标，经纬度；
 xk，yk          当前潜源顶点的坐标，以场点为原点，以当前破裂方向为x轴正方向，以km为单位
 wx,wy           潜源顶点重新排列后的坐标，以场点为原点，断层破裂方向为x轴正方向，断层距最小的顶点排第一位
 nang            每个潜源的顶点个数
 hk(20)         ph(20)  hk长轴衰减方向，ph衰减方向的概率
 zz              沿场点所在纬度小圆，移动1度的实际距离(km)
%}
%     ----------------------------------------------------------------
%% 定义全局变量 setup global variables
global   C CR R 

C=0.01745329252;  % c=pi/180;
R=6371.0;   % 地球平均半径
CR=C*R;     %  cr 是沿经线方向（南北方向）变化1度的距离km    


tic;
% background dirtribution ? yes;byes=1; 
nrd=fopen(inputfile,'rt');
% riskd=[0.1, 0.05,  0.02, 0.01,0.005,0.002,0.001,5e-3,2e-3,1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5,1e-5,5e-6,2e-6,1e-6];
% riskd=riskd';
%% 读入数据  import data from inputfile
% ============Read in GPV general parameter value  读入控制性参数=====================
[nsou, nsite, latt, modx,gsl,ca,cb,ti,riskd,sitex,sitey,sx,sy, dx, dy, nx, ny]=rdgpv(nrd);
% ca为长轴衰减系数；cb为短轴衰减系数；第一行为东部，第二行为西部

% ============Read in SPD seismic province data 读入地震带数据=====================
 [nz, nz1, amu, rate,beta,deps,amstep]=rdspd(nrd);
 
 
%=============Read in BRS background source data 读入背景数据==============================
background=rdbrs(nrd);


%=============Read in seismic source data 读入潜源数据======================================
[nq0,amq,ammin,af0,ae,nang,ammas,lst,hks,phs,xs,ys,x,y,q0]=rdsou(nrd,ti,nsou);

%% 计算各震级档的单位面积发生频率 Compute annual rate on unit area for magnitude intervals 

[ff]=vm(beta,amu,ammin,nsou,nz,nz1,amstep,rate,nq0,amq,q0,ae);


%% 计算场点个数 compute number of sites
% nx，ny表示经线方向和纬线方向的格点数
 if gsl==1         % 场点按网格排列方式输入       
     if nx==0
         nx=1;
     end
     if ny==0
        ny=1;
     end
     nsite=nx*ny;
 end
%========场点坐标及需要计算的潜源=======================================

%% 计算每个场点的地震危险性曲线  calculate seismic hazard for every site
for is=1:nsite
    if gsl==2    % 逐个输入?
        sstx=sitex(is);
        ssty=sitey(is);
    else
        js=floor((is-1)/nx);  %此处用场点数is减去1，后除以每行格点数，
        % 然后取其整数部分，得出场点所在行，计算ssty时，不需要减
        ks=is-js*nx;
        sstx=sx+(ks-1)*dx;
        ssty=sy+js*dy;
    end
     zz=cos(C*ssty)*CR;
    % 确定对场点is有影响的潜源 即AR (ammin, maxd)>=af0 
    % 确定场点位置，东西部
    if sstx >= 105
            iat=1;
    else
            iat=2;
    end
   [nl, ord]=involve(nsou,ca,ammas,iat,af0,CR,zz,xs,ys,sstx,ssty);
    % 若所有潜源对场点均无贡献，则跳出当前循环
    % 执行下一次循环，计算下个场点
    if nl==0; continue; end  
    [risk_median]=hazme(nz,nz1,nl,ord,nsou,sstx,ssty,nang,xs,ys,lst,ammas,ammin,ff,phs,hks,ca,cb,deps,ti,zz,amstep,modx);
	% trl : 衰减关系截断水平
    [risk]=correct(trl, risk_median,nl,ord,sigma,'normal',nz,nz1,nsou,sstx,ssty,nang,xs,ys,lst,ammas,ammin,ff,phs,hks,ca,cb,deps,ti,zz,amstep,modx);

    % 背景源的贡献, brs_cor为40*n_steps矩阵，每列为响应epsilon
    % 水平下背景源的贡献，计算中使用长轴衰减关系
    brs_cor=brs_correct(background,amstep,ca,cb, modx,ti, sigma,trl,'normal');
    brs=sum(brs_cor,2);
    risk=risk+byes*brs;
    att=exp(ti); 
    att=att';
    risk=1-exp(-risk);
    
    
	y=log(risk(risk>0));
    x=log(att(risk>0)); 
    y1=log(riskd);
    a=max(y);
    b=min(y);
	start=find(y1<=a,1,'first');
	finish=find(y1>=b,1,'last');
	y1=y1(start:finish);
	x1=interp1(y,x,y1);
	gmpa=exp(x1);
    riskd=exp(y1);
	% save  sha_dis.mat   risk   risk_median   att riskd gmpa
end
loglog(att, risk)
title('Seismic hazard curve')
xlabel('PGA(gal)');
ylabel('APE');
grid on
fclose('all');
toc;



function [nsou, nsite, latt, modx,gsl,ca,cb,ti,riskd,sitex,sitey,sx,sy, dx, dy, nx, ny]=rdgpv(nrd)
%{
RDGPV: read cpsha input file and write output to file
transmite data to variable in MATLAB workspace 
%}
%open input & output files

% read in general parameter value
gpv=fscanf(nrd,'%g',[6 1]);
nsou=gpv(1);       % 潜源个数
nsite=gpv(2);      % 场点个数
latt=gpv(3);       % 预定地震动参数： 1为对数,2为正常值
modx=gpv(4);       % 衰减关系模型，pga或sa为2；断层源时为3    
gsl=gpv(5);        % 场点输入方式,1时逐个输入，2时网格输入      
% slll=gpv(6);  
if nsou>2100
        error('潜源总数大于2100,计算终止!');
end
      
if nsite>3000
        error('场点总数大于3000,计算终止!');
end
 
      ca=fscanf(nrd,'%g', [1  14]); 
	  ca=reshape(ca, [7,2])'; % ca  东部衰减关系 or 衰减关系
      cb=fscanf(nrd,'%g', [1  14]);
	  cb=reshape(cb, [7, 2])'; % cb  西部衰减关系 or 震级-破裂长度关系
      

      
      % ca长轴衰减关系，第一行为东部，第二行为西部
      % cb短轴衰减关系，第一行为东部，第二行为西部
	  temp1=ca; temp2=cb;
      ca=[temp1(1,:); temp2(1,:)]; cb=[temp1(2,:); temp2(2,:)];
      % save  ar_ce.mat   ca  cb
      lnatt=fscanf(nrd,'%g',[1 40]);
      if latt~=1    % 为预定地震动幅值取自然对数         
          ti=log(lnatt);
      else
          ti=lnatt;
      end
      
      
      tnrd=fscanf(nrd, '%g', [1 10]);
      riskd=tnrd';

    
if gsl==2
			siteloc=fscanf(nrd, '%g', [2 nsite]);
            % 注意，fscanf 按行读入数据，fprintf 按列输出数据
            % 因此，需要把数据矩阵转置
            sitex=siteloc(1,:); sitey=siteloc(2,:);
			sx=0; sy=0; nx=0; ny=0; dx=0; dy=0;
else
			 siteloc=fscanf(nrd,'%g',[6 1]);
             sx=siteloc(1); sy=siteloc(2); dx=siteloc(3); dy=siteloc(4); 
             nx=siteloc(5); ny=siteloc(6); sitex=0; sitey=0;
end

function [nz, nz1, amu, rate,beta,deps,amstep]=rdspd(nrd)
%RDSPD read in seismci province  data
      nz=fscanf(nrd, '%d', [1 1]);

      if nz>30
        error('地震带个数大于30,计算终止!');
      end
      nz1=fscanf(nrd, '%g', [nz 1]);
	  nzdata=fscanf(nrd,'%g',[4 nz]);
      amu=nzdata(1,:); 
	  rate=nzdata(3,:);
	  beta=nzdata(2,:); 
	  deps=nzdata(4,:);
      amstep=fscanf(nrd,'%g', [1 1]);


function [nq0,amq,ammin,af0,ae,nang,ammas,lst,hks,phs,xs,ys,x,y,q0]=rdsou(nrd,ti,nsou)
global  CR   
      nq0=fscanf(nrd,'%g',[1 1]);   % nq0 震级分档个数
      if nq0>10
        error('震级分档数大10,计算终止!');
      end
      amq=fscanf(nrd,'%g',[nq0,1]);  % amq为震级档情况，每个震级档的起始震级      
      amq(nq0+1)=amq(nq0)+9.0;   
      ammin=amq(1);
      af0=ti(1);                     % af0幅值初始化，即最小的预定幅值
      ae=zeros(nsou,1); nang=zeros(nsou,1); ammas=ae; lst=ae;
      hks=zeros(2,nsou); phs=hks; xs=zeros(10, nsou); ys=zeros(10, nsou);
      q0=zeros(10, nsou);
for jf=1:nsou
    sjf=fscanf(nrd,'%g', [3 1]);
    nang(jf)=sjf(1); ammas(jf)=sjf(2); lst(jf)=sjf(3);
    nn5=nang(jf);  % 当前潜源顶点个数，ammas(jf)震级上限
    if nn5>10 
        error('潜源节点总数大于10,计算终止!');
    end
    tla=fscanf(nrd,'%g',[1 4]);   % tla, temp variable of long axis: 
    % 长轴方向及其概率
    hks(1,jf)=tla(1)*pi/180; hks(2,jf)=tla(3)*pi/180;
    phs(:,jf)=[tla(2); tla(4)];
    q0(1:nq0,jf)=fscanf(nrd, '%g', [nq0, 1]);   % q0空间分布函数，各震级档地震带内归1
    locate=fscanf(nrd,'%g', [2 nn5]); 
    xs(1:nn5,jf)=locate(1,:); ys(1:nn5,jf)=locate(2,:);
    x=zeros(nn5,1); y=zeros(nn5,1);
    for kj=1:nn5
        x(kj)=CR*cos(ys(1,jf)*pi/180)*(xs(kj,jf)-xs(1,jf));
        y(kj)=CR*(ys(kj,jf)-ys(1,jf));
%     cr 是沿/纬度变化1度的距离，unit：km
%     x(kj) y(kj)为第kj号顶点相对于1号顶点的位置，或者说是以1号顶点为原点， 以其纬线和经线为横纵轴的直角坐标系内的坐标
%     x表示横坐标，即沿纬线方向的距离，因此为小圆，乘以纬度的cos(ys)     
%     若为断层源，即潜源为直线时，nn5=2，ae为断层线的长度   
    end
    
    if nn5==2
       ae(jf)=sqrt(x(2)*x(2)+y(2)*y(2));
    else
        for kj=2:nn5-1
            ae(jf)=ae(jf)+abs(x(kj)*y(kj+1)-x(kj+1)*y(kj));
        end
    end
   
    ae(jf)=0.5*ae(jf);
end 

      
      
      
 

 
function [ff]=vm(beta,amu,ammin,nsou,nz,nz1,amstep,rate,nq0,amq,q0,ae)
%VM compute annual frequency of earthquakes in seismic source zone jf,
%magnitude interval i
ff=zeros(nsou,25);  
% ff为潜源i,震级档j,单位面积上的地震年平均发生率
% 第一层为ji=1：nz按地震带编号循环
% 第二层为i=nb1:nz1(ji),即当前地震带内的潜在震源区编号
% 第三层为当前潜源下每个震级档
for ji=1:nz
   beta(ji)=beta(ji)*log(10);
   if ji==1
        nb1=1;
   else  
        nb1=nz1(ji-1)+1;
   end
   amx=amu(ji); 
   ami=ammin;
   pp1=-beta(ji)*(amx-ami);
   ak=1.0/(1.0-exp(pp1));
   for i=nb1:nz1(ji)  % nb1 是当前地震带第一个潜源的编号，nz1(ji) 是第ji个地震带内的潜源编号的最大值       
       nm=(amu(ji)-ammin)/amstep+0.01;
       nm=floor(nm);
       
       for j=1:nm
          
           amag=ami+j*amstep;
           
           amp=amag+amstep/2.0;  % amp：amplitude plus震级档j内的震级上限
           amm=amag-amstep/2.0;  % amm：amplitude minus震级档j内的震级下限
           
           pp2=-beta(ji)*(amm-ami);
           pp3=-beta(ji)*(amp-ami);
           ep2=exp(pp2);
           ep3=exp(pp3);
           pmag=ak*rate(ji)*(ep2-ep3);  
           % pmag 是当前地震带震级档（amagm-dm/2, amag+dm/2)的地震的年平均发生率,
           % rate(ji)是整个地震带地震年平均发生率,还需要与空间分布函数q0相乘
           for ki=1:nq0
               if amag<amq(ki+1) && amag>=amq(ki)
                   ff(i,j)=q0(ki,i)/ae(i);
               end
           end
           ff(i,j)=ff(i,j)*pmag; 
        end 
% q0(ki,i) 是第ki个震级档、第i个潜源的空间分布函数，ff(i,j)是潜源i内震级档j单位面积上的地震年平均发生率
% if 语句给出amag的地震带ji、潜源i、震级档j，在整个震级档矩阵中的位置 知道了这个位置ki，才能引用空间分布函数q0(ki,i)
% 再用空间分布函数除以潜在震源区的面积


       
   end
end
%============================单位面积上的地震年平均发生率===================  
 
function  [l, ord]=involve(nsou,ca,ammas,iat,af0,CR,zz,xs,ys,sstx,ssty)
%INVOLE: PARTICATE IN COMPUTE FOR SITE IS
% MATLAB主要是传值，因此需要注意参数的顺序
    l=0;ord=zeros(nsou,1);
      ca1=ca(:,1);
      ca21=ca(:,2);
      ca22=ca(:,3);
      ca3=ca(:,4);
      ca8=ca(:,5);
      ca9=ca(:,6);
 
 for i=1:nsou
               
        ca2=ca21(iat)+ca22(iat)*ammas(i);   % ammas(i)是第i个潜源震级上限        
        ra0=ca8(iat)*exp(ca9(iat)*ammas(i));
        % ra是衰减所能及的最大距离        
        ra=exp((ca1(iat)+ca2*(ammas(i)-0.2)-af0)/ca3(iat))-ra0;
        % 转换为经纬度
        az1=ra/CR;  % az1是沿经线方向衰减所及的最大距离(用纬度表示）
        az2=ra/zz;  % az2是沿纬线方向衰减所及的最大距离(用经度表示）
        if az1>5;   az1=5.0; end
        if az2>5;   az2=5.0; end
        
        % 确定潜源经度范围x2-x1,纬度范围y2-y1
        x1=min(xs(:,i)); x2=max(xs(:,i));
        y1=min(ys(:,i)); y2=max(ys(:,i));
      
        % 确定场点is是否受当前潜源i的影响       
        if (abs(x1-sstx)<az2 || abs(x2-sstx)<az2) && (abs(y1-ssty)<az1 || abs(y2-ssty)<az1)
            % 确定潜源编号
            l=l+1;
            ord(l)=i;
        end
 end
 
 
 function  background=rdbrs(nrd)
 %RDBRS read background seismic source
       background=fscanf(nrd,'%g', [5 1]);

 function   brs=brs_dis(background,amstep,ca,cb, modx,ti)  
 %==============================背景源贡献================================
 % compute annual probability for excedence distribute from back ground
 % sources
 % 与场点位置无关，只与地震动幅值有关，影响场椭圆面积*单位面积年平均发生率
      bmu=background(1);
      bbeta=background(2); 
	  brate=background(3);
      bm0=background(4);
	  bdep=background(5);
      bbeta=-log(10)*bbeta;
      ek=exp(bbeta*(bmu-bm0));  
      ek=(1-ek)^(-1);  
%  ek  计算地震发生频率时的分母，即根据g-r关系，积分时使概率为1的常数

% ca和cb分别为短轴和长轴的衰减关系常数
      ca1=ca(:,1);
      ca21=ca(:,2);
      ca22=ca(:,3);
      ca3=ca(:,4);
      ca8=ca(:,5);
      ca9=ca(:,6);
 % 下为短轴常数
      cb1=cb(:,1);
      cb21=cb(:,2);
      cb22=cb(:,3);
      cb3=cb(:,4);
      cb8=cb(:,5);
      cb9=cb(:,6);
      
      
      brs=zeros(40,1);     %  brs背景地震对场地的影响 
      iat=1;               %  此处iat指定衰减关系轴，此处为长轴
    for j=1:40
        nm=(bmu-bm0)/amstep;
		nm=floor(nm);   
        for jj=1:nm
% 计算jj震级档中间震级对场地在场地处引起地震动幅值为a的震中距,衰减关系取中值  
            amag=bm0+(jj-0.5)*amstep;       %  amag为实际计算过程中震级档jj的中间震级           
            ca2=ca21(iat)+ca22(iat)*amag;  
            cb2=cb21(iat)+cb22(iat)*amag;
            ra0=ca8(iat)*exp(ca9(iat)*amag);  % ra与rb为衰减关系中的r0
            rb0=cb8(iat)*exp(cb9(iat)*amag);  
            if  modx~=3
                ra=(ca1(iat)+ca2*amag-ti(j))/ca3(iat);
%   ti(j)是预设幅值矩阵，给出了预先给定的值，取对数后的结果，即衰减关系中的lna 
                ra=exp(ra)-ra0;
                rb=(cb1(iat)+cb2*amag-ti(j))/cb3(iat); 
                %  (cb1(iat)+cb2*amag)衰减关系的前3项之和               
                rb=exp(rb)-rb0;
                if (rb<bdep || ra<bdep)
                    s=0.0;          % 背景源影响场椭圆面积
                else
                    ra=sqrt(ra*ra-bdep*bdep);    %  ra，rb为震中距
                    rb=sqrt(rb*rb-bdep*bdep);
                    s=ra*rb*pi;            % s 背景源影响场面积
                end
            else
                rll=10^(cb1(iat)+cb2*amag); % rll中间震级对应的断层破裂长度                
                ra=(ca1(iat)+ca2*amag-ti(j))/ca3(iat);
                ra=exp(ra)-ra0;
                if ra<=bdep
                     s=0.0;
                else
                    ra=sqrt(ra*ra-bdep*bdep);
                    s=rll*2.*ra+pi*ra*ra;  % 断层源影响面为长方形加圆形，ra为震中距
                end
            end
            pm=bbeta*(amag-bm0);
            pm=exp(pm)*amstep*ek;
            brs(j)=brs(j)+s*brate*pm;
        end     
    end 
  
 %==============================背景源贡献===============================
function   brs_cor=brs_correct(background,amstep,ca,cb, modx,ti, sigma,trl, dt)
cl=ca(:,1); cs=cb(:,1);
% 地震动衰减关系校正，trl截断水平
e_min=-trl; 
e_max=trl;
e_step=(e_max-e_min)/20;
cm=e_min-0.5*e_step:e_step:e_max-0.5*e_step;    % cm, c_minus,衰减关系校正中的区间下限
% cc=e_min:e_step:e_max;                                        % cc，c_center, …… 中间值
cp=e_min+0.5*e_step:e_step:e_max+0.5*e_step;    % cp，c_plus,   …… 上限
cpv=cdf(dt,cp,0,1)-cdf(dt,cm,0,1);                            % 区间内正态分布累计值
cpv=cpv/(cdf(dt,e_max,0,1)-cdf(dt,e_min,0,1));         % 相对累计值
n_steps=(e_max-e_min)/e_step+1; 
n_steps=floor(n_steps);

% 校正后的背景源贡献
brs_cor=zeros(40, n_steps);

for kk=1:n_steps
     ca(:,1)=cl+(e_min+kk*e_step)*sigma; 
     cb(:,1)=cs+(e_min+kk*e_step)*sigma;
     brs_cor(:,kk)=brs_dis(background,amstep,ca,cb, modx,ti)*cpv(kk);
end 
 
 
 
 function    [risk_median]=hazme(nz,nz1,nl,ord,nsou,sstx,ssty,nang,xs,ys,lst,ammas,ammin,ff,phs,hks,ca,cb,deps,ti,zz,amstep,modx)
%HAZME compute median value of seismic hazard
%{
% nl对场点is有贡献的潜源个数
% ord对场点有贡献的潜源编号数组，实际长度为nl
% ord2对场点有贡献的潜源应用衰减关系标度，实际长度为nl
% zz 是沿场点纬线（东西方向），经度变化一度的距离（km) 
%}


global CR C 
   
   ca1=ca(:,1); ca21=ca(:,2); ca22=ca(:,3); ca3=ca(:,4); ca8=ca(:,5); ca9=ca(:,6);
   cb1=cb(:,1); cb21=cb(:,2); cb22=cb(:,3); cb3=cb(:,4); cb8=cb(:,5); cb9=cb(:,6);

    ord2=zeros(nsou,1); 
    % 当前场点is的地震危险性矩阵, 40个地震动水平，每个cell元素为n_震级档*n_潜源矩阵
	risk_median=zeros(40, 1);

    %  默认考虑背景源
        
	%5重循环计算出当前场点地震危险性矩阵
	for lf=1:nl    %第一层循环，潜源jf=ord(lf)
        % 计算当前场点的地震危险性
        %每个潜源在场点处的地震动影响超过给定幅值的概率	
        jf=ord(lf);
        % 找出第jf个潜源所在的地震带      
        for kk=1:nz
            %下述if语句中一定要用<=, 否则会引起程序误跳出
            if jf<=nz1(kk) ; dep=deps(kk); break; end
            % break 跳出当前循环，执行循环外的下一条命令
        end
        nn5=nang(jf); 
        % 当前潜源jf经纬度坐标 
        x=xs(1:nn5,jf); y=ys(1:nn5,jf);
        % xav当前潜源顶点经度平均值
		xav=mean(x);
		if xav>=105 
            iat=1;
		else
            iat=2;
		end
        % 潜源衰减关系
        ord2(lf)=iat;   % 第lf个对当前场点有贡献的潜源，应用的衰减关系
        % 潜源顶点相对于场点的坐标，单位为km
        x(nn5+1)=x(1); y(nn5+1)=y(1);
		for i=1:nn5+1
            x(i)=-(sstx-x(i))*zz;
            y(i)=-(ssty-y(i))*CR;
		end
		
		d=1; % 潜源平行于断层破裂方向展布的条带的宽度，计算中的步长
        if  lst(jf)~=3             % 根据潜源类型，确定断层破裂方向及其概率            
            nd=2; hk=zeros(nd,1); ph=hk;
            for mm=1:nd
                hk(mm)=hks(mm,jf);   % hk为断层破裂方向延顺时针方向与N的夹角，即pi/2-hk(mm)
                ph(mm)=phs(mm,jf);
            end
        else
            nd=15; hk=zeros(nd,1); ph=hk;
			for mm=1:15
                hk(mm)=(mm-0.5)*12.0*C;
                ph(mm)=1./15.0;
			end
        end
		
        % mm=1:nd, 断层破裂方向，每个方向算一次
        % 乘以相应的概率
		for  mm=1:nd  % 第二层循环，
			if ph(mm)<=0 ; continue ; end
			ck1=sin(hk(mm));
			ck2=cos(hk(mm));
			% 坐标轴旋转，以断层破裂方向为x轴正方向，场点为原点
			xk=zeros(nn5+1,1); yk=xk;
			for i=1:nn5+1
				xk(i)=x(i)*ck2+y(i)*ck1;
				yk(i)=y(i)*ck2-x(i)*ck1;
			end
			if nn5==2 
				x1=xk(1); x2=xk(2); yb=abs(yk(1)); lk=1;
			end   
			% lk 潜源条带个数,平行于断层破裂方向			
			[no, lmin]=min(yk); [no1, l1]=max(yk); lk=floor((yk(l1)-yk(lmin))/d); 
			% 重新排列潜源顶点坐标，y值最小的排第1
			wx=zeros(nn5+1,1); wy=wx;
			for i=lmin:lmin+nn5
				j1=i-lmin+1;
				if i>nn5; j2=i-nn5; else j2=i; end
				wx(j1)=xk(j2); wy(j1)=yk(j2);
			end
			wx(nn5+1)=wx(1); wy(nn5+1)=wy(1);
			% g相邻顶点所确定直线的正切值，b,y轴截距	
			g=zeros(nn5,1); b=g;
			for i=1:nn5
				if  abs(wx(i+1)-wx(i))<0.005
					g(i)=500;
				else
					g(i)=(wy(i+1)-wy(i))/(wx(i+1)-wx(i));
				end
				b(i)=wy(i)-g(i)*wx(i);
			end
			for ji=1:lk     %第三层循环，计算每个条带对场点的影响
				y1=wy(1)+(ji-0.5)*d;  % y1场点垂直条带距离平均值
                for ks=nn5:-1:2
                    if y1<wy(ks) 
						b1=b(ks); g1=g(ks); break;
                    end
                end
				
                for ks=2:nn5
                    if y1<wy(ks) 
						b2=b(ks-1); g2=g(ks-1); break;
                    end
                end
				
                % x1,x2是直线y=y1与潜源边界交点的横坐标     
				x1=(y1-b1)/g1;   
				x2=(y1-b2)/g2;
				ammax=ammas(jf);
                % 当前潜源震级档个数，可能不是整数
				nmag=(ammax-ammin+0.1)/amstep;   
				nmag=floor(nmag+0.01);
				fsl=abs(x2-x1);
				yb=abs(y1);         % 约为场点到潜源条带最小距离的平均值
				xl1=abs(x1);
				xl2=abs(x2);
              
				for  i=1:40              % 第四层循环，每个地震动强度算一次
                 
					% 最小的ra0和rb0
					rb0=cb8(iat)*exp(cb9(iat)*ammin); 
					ra0=ca8(iat)*exp(ca9(iat)*ammin);
					if modx ~=3
						ammi=(ti(i)+cb3(iat)*log(yb+rb0)-cb1(iat))/cb21(iat)/1.2; 
					else
						ammi=(ti(i)+ca3(iat)*log(yb+ra0)-ca1(iat))/ca21(iat)/1.2;
					end
					kam=floor((ammi-ammin+0.1)/amstep);
					if kam<1; kam=1; end
					rpp=zeros(nmag,1);
					for im=kam:nmag     % 第五层循环，震级档
						rpp(im)=d*ff(jf,im)*ph(mm); 
						a=floor(im); 
						amag=ammin+(a-0.4)*amstep;
                        if modx~=3
                            ca2=ca21(iat)+ca22(iat)*amag;
                            cb2=cb21(iat)+cb22(iat)*amag;
                            ra0=ca8(iat)*exp(ca9(iat)*amag);
                            rb0=cb8(iat)*exp(cb9(iat)*amag);
                            ca4=(ca1(iat)+ca2*amag-ti(i))/ca3(iat);
                            cb4=(cb1(iat)+cb2*amag-ti(i))/cb3(iat);
                            ra=exp(ca4)-ra0;
                            if ra<dep; continue; end
                            ra=sqrt(abs(ra*ra-dep*dep));
                            rb=exp(cb4)-rb0;
                            if rb<dep; continue; end
                            rb=sqrt(abs(rb*rb-dep*dep));
                            if yb>=rb; continue; end
                            xa=ra*sqrt(abs(rb*rb-yb*yb))/rb;  
                            % 用rb来表示条带上所有点的纵坐标
                            % xa为对场点影响不小于ti(i)的最大值
                            xl3=xl1+xl2-fsl;
                            xl3=abs(xl3);  
                            if  xl3<0.1   % x1与x2异号
                                xa1=min(xa,xl2);
                                xa2=min(xa,xl1);
                                xa=xa1+xa2;
                            else
                                xl0=max(xl1,xl2);
                                xa=min(xa,xl0);
                                xa=xa-xl0+fsl;
                            end
                        else          %圆形衰减，配合震级-破裂关系
                            ca4=(ca1(iat)+ca2*amag-ti(i))/ca3(iat);
                            ra0=ca8(iat)*exp(ca9(iat)*amag);
                            ra=exp(ca4)-ra0;
                            if ra<=dep; continue; end
                            ra=sqrt(abs(ra*ra-dep*dep));
                            if ra<yb; continue; end
                            ra=sqrt(abs(ra^2-yb^2));
                            rll=10^(cb1(iat)+cb2*amag);
                            xa=0.5*rll+ra;  % 双侧破裂
                        end
                        if xa<=0; continue; end
                        % 潜源ord(lf)=jf，震级档im对场点is关于地震动幅值exp(ti(i))的超越概率
  						risk_median(i)=risk_median(i)+rpp(im)*xa;  
					end
				end
			end
		end
     
	end  %上述5重循环计算出当前场点地震危险性矩阵的初始值
    % 用于下述校正过程的起始输入
    % 按有贡献的潜源逐个校正
    
function [risk]=correct(trl, risk_median,nl,ord,sigma,dt,nz,nz1,nsou,sstx,ssty,nang,xs,ys,lst,ammas,ammin,ff,phs,hks,ca,cb,deps,ti,zz,amstep,modx)
%CORRECT modify site seismic hazard matrix risk
% can use different distribution: dt, string variable
% dt：distribution term 分布参数, 如 'normal'
[nr,nc]=size(risk_median);

cl=ca(:,1); cs=cb(:,1);

% 地震动衰减关系校正系数cpv，correct probability value
% trl截断水平
e_min=-trl; e_max=trl; e_step=(e_max-e_min)/20;
cm=e_min-0.5*e_step:e_step:e_max-0.5*e_step;
cc=e_min:e_step:e_max;
cp=e_min+0.5*e_step:e_step:e_max+0.5*e_step;
cpv=cdf(dt,cp,0,1)-cdf(dt,cm,0,1);  
cpv=cpv/(cdf(dt,e_max,0,1)-cdf(dt,e_min,0,1));
n_steps=(e_max-e_min)/e_step+1; n_steps=floor(n_steps);
risk_cor=zeros(nr,n_steps);



for kk=1:n_steps
     ca(:,1)=cl+(e_min+kk*e_step)*sigma; 
     cb(:,1)=cs+(e_min+kk*e_step)*sigma;
     risk_epsilon=hazme(nz,nz1,nl,ord,nsou,sstx,ssty,nang,xs,ys,lst,ammas,ammin,ff,phs,hks,ca,cb,deps,ti,zz,amstep,modx);
	 risk_cor(:,kk)=risk_cor(:,kk)+risk_epsilon*cpv(kk);	
end
risk=sum(risk_cor,2);  % sum the risk_cor in rows, to give a column vector
