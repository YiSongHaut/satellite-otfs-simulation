function ch = leo_otfs_channel(cfg,pulseOverride) % 生成目标论文式(5)-(9)的路径模型及导频相关等效DD信道。
% 用法：ch=leo_otfs_channel('NTN-TDL-A','rectangular')按名称调用；也支持leo_otfs_channel(cfg)输入完整配置。
% MATLAB R2023b；不需要5G/通信/统计工具箱；所有索引字段均为从0开始的物理索引。
% 依据：目标论文 DOI 10.1109/TVT.2026.3730193；TR38.811 V15.3.0 §6.9.2。
% 重要：本文未公开DS、散射角、公共频偏补偿量和Nhat；本代码将这些作为显式配置，不能承诺复原图3数值。
% 此处是标准PDP约束的有限路径实现，不等同于完整nrTDLChannel的Jakes时域随机过程。
% 公共大尺度幅度吸收到接收功率/SNR标定；不另外模拟目标论文未给出分布参数的降雨和闪烁。
% 标准来源：https://atisorg.s3.amazonaws.com/archive/3gpp-documents/Rel15/ATIS.3GPP.38.811.V1530.pdf 。
% 默认卫星公共频偏已补偿；否则20GHz下公共频偏远超16点多普勒窗，不能直接塞入小保护区。
% ch.hDD为完整M乘N网格，ch.hWindow为原文截断窗口，ch.hFullDD为有物理路径的时延行乘全部N列；不逐次归一化。
if nargin<1 || ischar(cfg) || isstring(cfg) % 没有完整配置时允许直接按PDP名称调用单独的信道子函数。
    if nargin<1, selectedProfile='NTN-TDL-B'; else, selectedProfile=char(cfg); end % 默认B，也允许用户直接指定A/B/C/D。
    cfg=struct('M',32,'N',16,'deltaF',60e3,'fc',20e9,'satSpeed',7110,'altitude',1200e3,'elevationDeg',50,'userSpeed',500/3.6); % 原论文表I公开的通信系统和运动参数。
    cfg.profile=selectedProfile; cfg.delaySpread=1e-6; cfg.raysPerTap=16; % 未公开：默认DS选1us，散射子射线数选16。
    cfg.compensateSatellite=true; cfg.residualCommonHz=0; cfg.Nhat=1; cfg.lp=16; cfg.pulse='ideal'; % 未公开：公共频偏已补偿、有限泄漏宽度和导频位置。
end % 完成独立信道生成的缺省配置。
if nargin>=2, cfg.pulse=char(pulseOverride); end % 第二输入仅替换脉冲类型，保持同一组PDP和系统参数。
required = {'M','N','deltaF','fc','satSpeed','altitude','elevationDeg','userSpeed','delaySpread','profile','raysPerTap','compensateSatellite','residualCommonHz','Nhat','lp','pulse'}; % 必需配置名。
for j = 1:numel(required) % 检查接口，避免隐式使用未定义参数。
    assert(isfield(cfg,required{j}),['缺少cfg.',required{j}]); % 报告具体缺失字段。
end % 完成配置存在性检查。
assert(cfg.M>=2 && cfg.N>=2 && cfg.deltaF>0 && cfg.delaySpread>0,'网格或时延扩展非法'); % 检查基本物理尺度。
assert(cfg.raysPerTap>=1 && cfg.raysPerTap==round(cfg.raysPerTap),'raysPerTap须为正整数'); % 每个散射tap的子路径数。
assert(cfg.Nhat>=0 && cfg.Nhat==round(cfg.Nhat),'Nhat须为非负整数'); % 泄漏截断半宽。
assert(cfg.lp>=0 && cfg.lp<cfg.M,'lp必须是合法的0基延迟索引'); % 检查导频延迟坐标。
c0 = 299792458; % 光速，单位m/s。
earthRadius = 6371000; % 地球半径，单位m。
T = 1/cfg.deltaF; % 与论文理想临界采样约定一致，不额外加入CP长度。
switch upper(char(cfg.profile)) % 按名称调用独立的标准PDP子函数，方便逐个核对和扩展。
    case 'NTN-TDL-A' % 3GPP TR38.811表6.9.2-1，非视距。
        [delays,powersDb,isLos]=profile_ntn_tdl_a(); % 返回A的归一化时延、功率和LOS标志。
    case 'NTN-TDL-B' % 3GPP TR38.811表6.9.2-2，非视距。
        [delays,powersDb,isLos]=profile_ntn_tdl_b(); % 返回B的归一化时延、功率和LOS标志。
    case 'NTN-TDL-C' % 保持原有C模型入口供额外检验。
        [delays,powersDb,isLos]=profile_ntn_tdl_c(); % 返回C的归一化时延、功率和LOS标志。
    case 'NTN-TDL-D' % 3GPP TR38.811表6.9.2-4，视距。
        [delays,powersDb,isLos]=profile_ntn_tdl_d(); % 返回D的归一化时延、功率和LOS标志。
    otherwise % 拒绝未知PDP，避免不知不觉换成别的模型。
        error('profile必须是NTN-TDL-A/B/C/D'); % 提供清晰配置错误。
end % 结束PDP选择。
powerLin = 10.^(powersDb/10); powerLin = powerLin/sum(powerLin); % 归一化平均总路径功率，不是逐次归一化信道。
delaySeconds = delays*cfg.delaySpread; % 标准归一化时延乘目标RMS时延尺度。
delayTap = round(delaySeconds*cfg.M*cfg.deltaF); % 依论文整数时延假设映射；量化误差另行输出。
lmax = max(delayTap); % 全部路径的最大整数时延索引。
assert(lmax<cfg.M,'整数时延越界，请减小DS或增大M'); % 禁止无提示环绕物理时延。
satelliteHz = cfg.satSpeed/c0*cfg.fc*earthRadius/(earthRadius+cfg.altitude)*cosd(cfg.elevationDeg); % TR38.811的公共卫星频移。
commonHz = satelliteHz*(~cfg.compensateSatellite)+cfg.residualCommonHz; % 补偿后保留的公共频移。
ueMaxHz = cfg.userSpeed/c0*cfg.fc; % UE最大散射多普勒，单位Hz。
kmax = floor((abs(commonHz)+ueMaxHz)*cfg.N*T+0.5); % 覆盖所有最近整数多普勒中心，而非向上取整连续频移。
knu = kmax+cfg.Nhat; % 原文保护区半宽k_nu=kmax+Nhat。
assert(2*knu+1<=cfg.N,'观测多普勒窗大于N；须先补偿公共频偏或调整网格'); % 避免同一周期抽头重复。
gain = complex(zeros(0,1)); ell = zeros(0,1); nuHz = zeros(0,1); % 初始化可变数量的物理路径。
for a = 1:numel(delays) % 按每个PDP分量生成路径。
    nr = cfg.raysPerTap; % 散射tap默认细分的子路径数。
    if isLos(a), nr = 1; end % LOS镜面项不拆成多个独立镜面项。
    for r = 1:nr % 生成当前tap的每条子路径。
        if isLos(a) % LOS采用固定模长和随机初相。
            g = sqrt(powerLin(a))*exp(1i*2*pi*rand); theta = 0; % 这里的LOS到达方位为显式实现约定。
        else % 散射路径采用复高斯增益。
            g = sqrt(powerLin(a)/nr/2)*(randn+1i*randn); theta = 2*pi*rand; % 子路径共享PDP功率，散射方位均匀。
        end % 完成一条路径增益的生成。
        gain(end+1,1) = g; ell(end+1,1) = delayTap(a); nuHz(end+1,1) = commonHz+ueMaxHz*cos(theta); % 保存实际路径而不提前合并同格路径。
    end % 完成当前tap的路径生成。
end % 完成所有物理路径生成。
nuBin = nuHz*cfg.N*T; k = floor(nuBin+0.5); kappa = nuBin-k; % 将连续归一化多普勒拆分为整数与分数。
onBoundary = abs(kappa+0.5)<1e-12; k(onBoundary)=k(onBoundary)-1; kappa(onBoundary)=0.5; % 实现论文(-1/2,1/2]端点约定。
dAxis = -knu:knu; hDD = complex(zeros(lmax+1,numel(dAxis))); % 延迟为行、多普勒为列，hDD(:)就是延迟快变顺序。
hFullDD = complex(zeros(lmax+1,cfg.N)); % 保存未截断的整周期泄漏，用于量化截断误差。
dFullAxis = -floor(cfg.N/2):(ceil(cfg.N/2)-1); % 定义整周期多普勒显示坐标；偶数N时为-N/2到N/2-1。
for i = 1:numel(gain) % 将每条路径的泄漏贡献累加到等效信道。
    if strcmpi(cfg.pulse,'ideal') % 论文式(6)/(9)的双正交波形相位。
        phase = exp(-1i*2*pi*ell(i)*nuBin(i)/(cfg.M*cfg.N)); % 相位与接收延迟无关。
    elseif strcmpi(cfg.pulse,'rectangular') % 论文矩形波形在导频观测延迟l=lp+ell_i处的相位。
        phase = exp(1i*2*pi*cfg.lp*nuBin(i)/(cfg.M*cfg.N)); % 这里定义的是导频专用等效信道，不滥用于整帧数据卷积。
    else % 拒绝未实现波形。
        error('pulse必须为ideal或rectangular'); % 提示可用配置。
    end % 完成波形相位计算。
    for q = -cfg.Nhat:cfg.Nhat % 论文式(5)/(8)的有限泄漏范围。
        d = k(i)-q; kernel = mean(exp(1i*2*pi*(0:cfg.N-1)*(q+kappa(i))/cfg.N)); % 几何级数等价式，避免κ=0时的0/0。
        hDD(ell(i)+1,d+knu+1) = hDD(ell(i)+1,d+knu+1)+gain(i)*phase*kernel; % 相同DD格的贡献相干相加。
    end % 完成有限泄漏叠加。
    for dc = 0:cfg.N-1 % 完整周期的每个离散多普勒格。
        kernel = mean(exp(1i*2*pi*(0:cfg.N-1)*(nuBin(i)-dc)/cfg.N)); % 完整Dirichlet核，不能与截断真值混为一谈。
        hFullDD(ell(i)+1,dc+1) = hFullDD(ell(i)+1,dc+1)+gain(i)*phase*kernel; % 计算全周期参考信道。
    end % 完成全周期泄漏叠加。
end % 完成等效信道生成。
hTruncFull = complex(zeros(size(hFullDD))); hTruncFull(:,mod(dAxis,cfg.N)+1)=hDD; % 将有限窗口嵌入完整多普勒周期。
fullGrid = complex(zeros(cfg.M,cfg.N)); fullGrid(1:lmax+1,:)=hFullDD; % 将物理非零时延区嵌入完整M乘N离散DD网格。
windowGrid = complex(zeros(cfg.M,cfg.N)); windowGrid(1:lmax+1,:)=hTruncFull; % 将原文有限泄漏近似嵌入同样大小的DD网格。
ch.hDD=fullGrid; ch.h=hDD(:); ch.hWindow=hDD; ch.hFullDD=hFullDD; % 保留原估计向量接口，同时提供完整网格及窗口矩阵。
ch.hGridFinite=windowGrid; ch.hFullVector=hFullDD(:); ch.dFullAxis=dFullAxis; % 分别输出截断网格、完整物理向量和居中的多普勒轴。
ch.hCompleteActive=hFullDD(:,mod(dFullAxis,cfg.N)+1); ch.hGridVector=fullGrid(:); % 提供居中排列的完整物理系数和完整网格向量供先验建模。
ch.gain=gain; ch.delayTap=ell; ch.nuHz=nuHz; ch.nuBin=nuBin; ch.k=k; ch.kappa=kappa; % 输出调试用物理路径参数。
ch.lmax=lmax; ch.kmax=kmax; ch.knu=knu; ch.dAxis=dAxis; ch.T=T; % 输出系统和索引元数据。
ch.satelliteHz=satelliteHz; ch.commonHz=commonHz; ch.profile=cfg.profile; % 保存公共频移和PDP标签。
ch.delayQuantizationError=delaySeconds-delayTap/(cfg.M*cfg.deltaF); % 量化误差，单位秒。
ch.normalizedDelays=delays; ch.delaySeconds=delaySeconds; ch.profilePowerDb=powersDb; % 输出标准PDP和配置DS对应的实际时延以便追踪每行来源。
ch.profileIsLos=isLos; ch.profileQuantizedDelayTap=delayTap; ch.profilePowerLinear=powerLin; % 直接导出抽头级数据，不必从已细分的子射线反推。
ch.profileKFactorDb=NaN; % 非视距模型没有第一抽头的莱斯K因子。
if any(isLos), ch.profileKFactorDb=10*log10(sum(powerLin(isLos))/sum(powerLin(~isLos & delays==delays(find(isLos,1))))); end % 视距模型按原始零时延镜面/散射功率计算第一抽头K因子。
ch.truncationNMSE=norm(hFullDD-hTruncFull,'fro')^2/max(norm(hFullDD,'fro')^2,realmin); % 真实完整信道相对有限模型的误差。
ch.expectedPathPower=sum(powerLin); ch.numberPaths=numel(gain); % 平均功率为1；单次能量一般不为1。
end % 结束独立信道生成函数。

function [delays,powersDb,isLos]=profile_ntn_tdl_a() % 非视距A：三个Rayleigh抽头。
delays=[0;1.0811;2.8416]; powersDb=[0;-4.675;-6.482]; isLos=false(3,1); % 3GPP标准表6.9.2-1的原始数据。
end % 结束NTN-TDL-A抽头配置。

function [delays,powersDb,isLos]=profile_ntn_tdl_b() % 非视距B：三个较早抽头和一个较晚的弱抽头。
delays=[0;0.7249;0.7410;5.7392]; powersDb=[0;-1.973;-4.332;-11.914]; isLos=false(4,1); % 3GPP标准表6.9.2-2的原始数据。
end % 结束NTN-TDL-B抽头配置。

function [delays,powersDb,isLos]=profile_ntn_tdl_c() % 视距C：第一抽头含镜面和散射分量。
delays=[0;0;14.8124]; powersDb=[-0.394;-10.618;-23.373]; isLos=[true;false;false]; % 3GPP标准表6.9.2-3的原始数据。
end % 结束NTN-TDL-C抽头配置。

function [delays,powersDb,isLos]=profile_ntn_tdl_d() % 视距D：第一抽头含镜面和散射分量，另有两条时延路径。
delays=[0;0;0.5596;7.3340]; powersDb=[-0.284;-11.991;-9.887;-16.771]; isLos=[true;false;false;false]; % 3GPP标准表6.9.2-4的原始数据。
end % 结束NTN-TDL-D抽头配置。
