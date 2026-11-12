function results = run_otfs_nmse(mode,overrides) % 运行原文DD导频观测下的NMSE-SNR实验，各算法为本文件局部function。
% 使用：把两个m文件放同一目录，运行r=run_otfs_nmse('smoke')先检查；正式运行r=run_otfs_nmse('full')。
% 自定义：r=run_otfs_nmse('full',struct('nMC',500,'profile','NTN-TDL-A'))。
% 六情景对照：s=run_otfs_nmse('compareSmoke')先检查；s=run_otfs_nmse('compare',struct('nMC',500))再做正式统计。
% 本文件仅复现导频辅助信道估计，不加入UAMP/BG，不实现BER和时域采样接收机。
% 严格公式来源：NNCP-VSBL按目标论文式(19)-(36)，默认循环一维邻接保持原文写法。
% NCP-VB-grid按Dai等作者公开Burst_CE.m的固定网格先验更新，不包括其角度网格细化。
% 来源：https://github.com/jsdaiustc/Burst_CE/blob/main/Burst_CE.m ，访问日期2026-10-30。
% VB-SBL是独立Gamma-Gaussian复数均值场基线，不冒称已逐式复现Wang2022完整联合检测算法。
% SBL-EM是独立高斯方差的标准EM基线；OMP是固定预算贪婪基线；LS用于模型和SNR自检。
% Fig.3的BSBL采用分块相关高斯先验及在线块重组；APCSBL采用邻接模式耦合先验及数据驱动的双侧耦合强度。
% 两个参考算法的完整更新式和参数未在目标论文中列出；当前实现为结构一致的EM复现，不宣称与原作者代码逐式一致。
% 论文未公开的DS、Nhat、导频功率、SNR定义、初始化均在下面标为“实现选择”；不能据此宣称重现图3绝对数值。
% paperWindow保留式(5)的有限泄漏与原论文比较；complete用完整16列泄漏+纯导频观测探索全信道先验，二者不是同一实验。
% MATLAB R2023b兼容；仅用基础MATLAB函数；每行代码附中文说明；输出MAT/PNG/FIG由用户运行时生成。
% 验证边界：已做MATLAB语法解析及独立Python数值核验；交付环境无MATLAB，尚未实机执行本文件。
if nargin<1, mode='smoke'; end % 缺省先用小样本检查流程，不能把smoke曲线作为论文结论。
if nargin<2, overrides=struct; end % 可通过结构体覆盖实验配置。
if strcmpi(mode,'compare') || strcmpi(mode,'compareSmoke') % 用同一入口调用A/B/D乘双正交/矩形共六种情景。
    if strcmpi(mode,'compareSmoke'), caseMode='smoke'; else, caseMode='full'; end % 快速模式仅用于查错，正式模式按用户指定的MC统计。
    results=compare_six_scenarios(caseMode,overrides); return; % 六情景的原始结果和总览由比较子函数返回。
end % 完成比较入口的分流。
cfg.M=32; cfg.N=16; cfg.deltaF=60e3; cfg.fc=20e9; % 目标论文表I：子载波数、时隙数、子载波间隔、载频。
cfg.satSpeed=7110; cfg.altitude=1200e3; cfg.elevationDeg=50; cfg.userSpeed=500/3.6; % 表I：卫星和用户速度为m/s。
cfg.profile='NTN-TDL-B'; cfg.delaySpread=1e-6; cfg.raysPerTap=16; % 原文未公开DS和散射射线数：1us减少时延量化塌缩，16条射线近似标准各tap的散射多普勒谱；均须做敏感性扫描。
cfg.compensateSatellite=true; cfg.residualCommonHz=0; % 实现选择：完成公共卫星频偏预补偿，残余公共频偏可调。
cfg.Nhat=1; cfg.pulse='ideal'; cfg.lp=16; cfg.kp=8; % 实现选择：有限泄漏半宽1、双正交脉冲和0基导频位置。
cfg.estimationExtent='paperWindow'; % 默认沿用原文有限泄漏估计；设为complete时改用整周期泄漏的纯导频实验。
cfg.pilotEnergy=1; cfg.snrDefinition='data'; % 实现选择：SNR=单位平均QPSK数据功率/噪声功率，不逐次按接收能量调噪声。
cfg.snrDb=-5:5:25; cfg.nMC=100; cfg.seed=20261030; % 正式扫描范围、独立信道实现数和可重复随机种子。
cfg.maxIter=2500; cfg.minIter=10; cfg.tol=1e-5; cfg.floor=1e-12; % 迭代上限按收敛诊断提高；每次运行同时报告实际收敛率。
cfg.fig3Algorithms=false; % Fig.3图例算法由比较器打开，单情景旧入口仍可保留原基线。
cfg.initialAlpha=1; cfg.initialGamma=1; cfg.initialTau=1; % 实现选择：论文未公开的初始均值，不能借用真实信道初始化。
cfg.neighborhood='paperCircular'; cfg.ompAtoms=4; % 原文一维循环邻接；OMP固定4个原子而非偷用真实稀疏度。
cfg.outputDir=fullfile(pwd,'otfs_nmse_results'); cfg.saveOutputs=true; cfg.makePlots=true; cfg.verbose=true; % 结果保存、绘图和逐次进度开关。
if strcmpi(mode,'smoke') % 小规模调试模式只验证运行流程。
    cfg.nMC=3; cfg.snrDb=[0,10,20]; cfg.maxIter=40; % 减少蒙特卡洛和迭代量，不代表正式统计结果。
elseif ~strcmpi(mode,'full') % 防止误拼模式名称。
    error('mode应为smoke或full'); % 给出可选运行模式。
end % 完成模式配置。
names=fieldnames(overrides); % 读取用户显式覆盖的配置。
for j=1:numel(names) % 对每一个覆盖字段应用参数。
    assert(isfield(cfg,names{j}),['未知配置字段：',names{j}]); cfg.(names{j})=overrides.(names{j}); % 拒绝拼写错误造成的静默忽略。
end % 完成用户配置覆盖。
assert(cfg.nMC>=1 && cfg.nMC==round(cfg.nMC),'nMC必须是正整数'); % 检查实验次数。
assert(cfg.kp>=0 && cfg.kp<cfg.N && cfg.pilotEnergy>0,'导频配置非法'); % 检查导频坐标和功率。
assert(cfg.lp>=0 && cfg.lp<cfg.M,'导频延迟坐标越界'); % 检查另一个导频坐标。
rng(cfg.seed,'twister'); % 固定随机数种子，所有算法和SNR共用配对信道实现。
if cfg.fig3Algorithms % Fig.3采用图例中的六个方案名称和顺序。
    labels={'OMP','BSBL','VSBL','APCSBL','NCP-VSBL','NNCP-VSBL'}; % Fig.3图例顺序。
    algorithms={@omp_est,@bsbl_br,@vb_sbl,@apcsbl_adaptive,@ncp_vb,@nncp_vsbl}; % 六种算法均可独立调用。
else % 保留既有脚本原来的基线入口。
    labels={'LS','OMP-fixedK','SBL-EM','VB-SBL','NCP-VB-grid','NNCP-VSBL'}; % 旧模式不会误标为论文Fig.3。
    algorithms={@ls_est,@omp_est,@sbl_em,@vb_sbl,@ncp_vb,@nncp_vsbl}; % 原有算法集合。
end % 完成算法集合选择。
nAlg=numel(labels); nSNR=numel(cfg.snrDb); % 提前确定结果数组维度。
nmse=zeros(cfg.nMC,nSNR,nAlg); runtime=nmse; iterations=nmse; finalChange=nan(size(nmse)); converged=false(size(nmse)); % 保存逐次误差、时间、最终相对变化和收敛情况。
truncation=zeros(cfg.nMC,1); modelError=zeros(cfg.nMC,1); % 单独保存截断误差与前向模型闭环误差。
delayPower=zeros(cfg.nMC,cfg.M); % 每次记录完整16列信道在各离散时延抽头上的总能量。
noiseVariance=zeros(cfg.nMC,nSNR); estimatedAlpha=nan(size(nmse)); estimatedTau=nan(size(nmse)); % 保存噪声和学习到的精度。
diagnostics=cell(nSNR,nAlg); % 保存首个信道在各SNR下的完整算法轨迹供检查。
opts=struct('maxIter',cfg.maxIter,'minIter',cfg.minIter,'tol',cfg.tol,'floor',cfg.floor,'alpha0',cfg.initialAlpha,'gamma0',cfg.initialGamma,'tau0',cfg.initialTau,'ompAtoms',cfg.ompAtoms); % 初始化共用迭代参数，不预设与实际信道维度不符的连续块索引。
fprintf('模式=%s；估计支持=%s；PDP=%s；DS=%.3g us；射线/tap=%d；脉冲=%s；MC=%d。\n',mode,cfg.estimationExtent,cfg.profile,cfg.delaySpread*1e6,cfg.raysPerTap,cfg.pulse,cfg.nMC); % 输出原文公开参数与实现选择。
for trial=1:cfg.nMC % 每次生成一个独立物理信道，并在所有SNR上复用。
    ch=leo_otfs_channel(cfg); % 调用另一个m文件，生成真实路径、有限泄漏与等效信道。
    isComplete=strcmpi(cfg.estimationExtent,'complete'); % 完整泄漏实验的观测和先验向量长度均改为完整物理非零时延区。
    assert(isComplete || strcmpi(cfg.estimationExtent,'paperWindow'),'estimationExtent须为paperWindow或complete'); % 禁止未知模型名称。
    if isComplete % 完整模型与有限模型分别定义真值，避免拿不同目标计算NMSE。
        estAxis=ch.dFullAxis; ch.h=reshape(ch.hFullDD(:,mod(estAxis,cfg.N)+1),[],1); % 按居中多普勒列重排，按列堆叠全部16个离散多普勒抽头。
    else % 原文有限泄漏的局部估计模式。
        estAxis=ch.dAxis; ch.h=ch.hWindow(:); % 保持原文L=(lmax+1)(2kν+1)以及原比较算法入口。
    end % 完成估计对象的定义。
    assert(isComplete || 4*ch.knu+1<=cfg.N,'原文单导频保护区装不下：请减Nhat或调整网格'); % 有限窗口仍须容纳双侧保护区。
    assert(2*ch.lmax+1<=cfg.M,'单导频时延保护区装不下：请减小DS或增大M'); % 保证物理时延不环绕导频。
    assert(cfg.lp>=ch.lmax && cfg.lp+ch.lmax<cfg.M,'须令lp>=lmax且lp+lmax<M以使用矩形脉冲CE表达式'); % 保证CE观测区的相位公式有效。
    if isComplete % 整周期分数多普勒泄漏可到达全部N列，原文局部保护区已不足以隔离数据。
        xDD=complex(zeros(cfg.M,cfg.N)); % 纯导频观测用于无偏研究完整等效信道；此模式不与原论文的导频开销直接比较。
    else % 保留论文式(5)的有限泄漏假设和数据传输。
        xDD=((2*randi([0,1],cfg.M,cfg.N)-1)+1i*(2*randi([0,1],cfg.M,cfg.N)-1))/sqrt(2); % 生成单位平均功率4QAM数据。
        guardRows=mod(cfg.lp+(-ch.lmax:ch.lmax),cfg.M)+1; guardCols=mod(cfg.kp+(-2*ch.knu:2*ch.knu),cfg.N)+1; % 图1的时延和多普勒保护区。
        xDD(guardRows,guardCols)=0; % 严格按局部窗口清空保护区而不把全帧数据删去。
    end % 完成对应模式的导频与数据设计。
    xDD(cfg.lp+1,cfg.kp+1)=sqrt(cfg.pilotEnergy); % 在0基位置(lp,kp)插入唯一导频。
    xTF=fft(ifft(xDD,[],2),[],1)*sqrt(cfg.N/cfg.M); % 论文式(1)的ISFFT，矩阵行为延迟/子载波，列为多普勒/时间。
    xBack=ifft(fft(xTF,[],2),[],1)*sqrt(cfg.M/cfg.N); % 对应式(4)的逆变换自检，不假装实现了时域Heisenberg接收机。
    assert(norm(xBack-xDD,'fro')<=1e-10*max(norm(xDD,'fro'),1),'ISFFT/SFFT尺度错误'); % 校验系统变换归一化。
    ySignal=apply_channel_dd(xDD,ch,cfg); % 严格从路径和式(5)生成接收信号，不直接用Phi*h替代信道传输。
    obsRows=cfg.lp+(0:ch.lmax)+1; obsCols=mod(cfg.kp+estAxis,cfg.N)+1; % 观测时延窗相同，多普勒窗依模式分别为局部或全部N列。
    yClean=reshape(ySignal(obsRows,obsCols),[],1); % 按延迟快变的共同顺序堆叠观测。
    A=build_measurement(xDD,ch,cfg,estAxis); % 按各自多普勒支持独立构建对应测量矩阵。
    modelError(trial)=norm(yClean-A*ch.h)/max(norm(yClean),cfg.floor); % 比较物理前向计算与估计问题模型。
    assert(modelError(trial)<1e-9,'yClean与Phi*h不一致：检查相位、索引和保护区'); % 模型不一致时停止，绝不带错跑曲线。
    opts.nDelay=ch.lmax+1; opts.nDoppler=numel(estAxis); % h按时延优先、每个多普勒列含nDelay个系数。
    opts.neighbors=make_neighbors(opts.nDelay,opts.nDoppler,cfg.neighborhood); % 可在论文展平邻接和同一时延内多普勒邻接之间切换。
    unitNoise=(randn(cfg.M,cfg.N)+1i*randn(cfg.M,cfg.N))/sqrt(2); % 一份单位方差复高斯噪声在各SNR缩放，全部算法共用。
    truncation(trial)=ch.truncationNMSE; % 报告有限Nhat相对完整泄漏的误差，NMSE真值仍为同一有限h。
    delayPower(trial,1:ch.lmax+1)=sum(abs(ch.hFullDD).^2,2).'; % 保存同一完整信道的延迟边际功率用于跨模型比较。
    for is=1:nSNR % 遍历SNR轴。
        if strcmpi(cfg.snrDefinition,'data') % 默认按数据功率定义SNR，保持Rayleigh深衰落。
            sigma2=10^(-cfg.snrDb(is)/10); % 数据平均符号功率为1。
        elseif strcmpi(cfg.snrDefinition,'observation') % 备选逐次观测平均功率SNR，必须与默认结果区分标签。
            sigma2=mean(abs(yClean).^2)*10^(-cfg.snrDb(is)/10); % 注意这一设置会消除部分深衰落影响。
        else % 阻止未知SNR定义。
            error('snrDefinition应为data或observation'); % 提示合法模式。
        end % 完成噪声功率确定。
        yNoisy=ySignal+sqrt(sigma2)*unitNoise; y=reshape(yNoisy(obsRows,obsCols),[],1); % 在物理接收网格加噪，再提取观测。
        noiseVariance(trial,is)=sigma2; % 记录真实噪声，仅用于事后检查而不提供给估计器。
        for ia=1:nAlg % 在完全相同的A和y上运行各算法。
            clockStart=tic; [hhat,info]=algorithms{ia}(y,A,opts); runtime(trial,is,ia)=toc(clockStart); % 调用独立function并计时。
            assert(all(isfinite(hhat)),'估计器出现NaN/Inf，停止并排查，不静默删除失败样本'); % 数值失败不做有利筛选。
            nmse(trial,is,ia)=norm(hhat-ch.h)^2/max(norm(ch.h)^2,cfg.floor); % 对每个独立信道计算相对误差，在线性域平均。
            iterations(trial,is,ia)=info.iterations; converged(trial,is,ia)=info.converged; if isfield(info,'finalChange'), finalChange(trial,is,ia)=info.finalChange; elseif isfield(info,'traceChange'), finalChange(trial,is,ia)=info.traceChange(end); end % 记录迭代数、收敛标记和最终相对变化。
            estimatedAlpha(trial,is,ia)=info.alpha; estimatedTau(trial,is,ia)=info.tau; % 保存精度学习诊断。
            if trial==1, diagnostics{is,ia}=info; end % 首实例保留phi、gamma和alpha/tau收敛轨迹。
        end % 完成当前SNR的所有算法。
    end % 完成当前信道的所有SNR。
    if trial==1 % 首次保存完整系统实例方便断点检查。
        example=struct('channel',ch,'xDD',xDD,'xTF',xTF,'A',A,'yClean',yClean,'obsRows',obsRows,'obsCols',obsCols,'estAxis',estAxis); % 保存完整信道与所用观测轴。
        example.gramCondition=cond(A'*A); example.scaledIdentityError=norm(A-sqrt(cfg.pilotEnergy)*eye(size(A)),'fro'); % 检查单导频字典是否退化为缩放单位阵。
    end % 完成首实例保存。
    if cfg.verbose, fprintf('已完成 %d/%d 个信道，模型误差 %.2e，泄漏截断NMSE %.3g。\n',trial,cfg.nMC,modelError(trial),truncation(trial)); end % 比较模式可关闭逐次打印。
end % 完成全部蒙特卡洛实验。
meanNMSE=nan(nSNR,nAlg); se=nan(nSNR,nAlg); validCount=zeros(nSNR,nAlg); % 统计时只将已收敛估计计入主NMSE曲线。
for is=1:nSNR, for ia=1:nAlg, valid=converged(:,is,ia); validCount(is,ia)=sum(valid); if any(valid), vals=nmse(valid,is,ia); meanNMSE(is,ia)=mean(vals); if numel(vals)>1, se(is,ia)=std(vals,0)/sqrt(numel(vals)); else, se(is,ia)=NaN; end; end; end; end % 未收敛样本保留在rawNMSE中，但不静默纳入主均值。
meanNMSEAll=reshape(mean(nmse,1),nSNR,nAlg); % 另存包含全部样本的原始均值以便审计。
meanNMSEDb=nan(size(meanNMSE)); hasValid=isfinite(meanNMSE); meanNMSEDb(hasValid)=10*log10(max(meanNMSE(hasValid),cfg.floor)); ciLow=nan(size(meanNMSE)); ciHigh=nan(size(meanNMSE)); ciLow(hasValid)=max(meanNMSE(hasValid)-1.96*se(hasValid),0); ciHigh(hasValid)=meanNMSE(hasValid)+1.96*se(hasValid); % 零个收敛样本时保留NaN，不把无结果伪装成地板NMSE。
results=struct('cfg',cfg,'labels',{labels},'nmseLinear',meanNMSE,'nmseDb',meanNMSEDb,'rawNMSE',nmse,'nmseLinearAll',meanNMSEAll); % 组织收敛样本主结果和含所有样本的审计结果。
results.validCount=validCount; results.convergedFraction=validCount/cfg.nMC; results.ciLinearLow=ciLow; results.ciLinearHigh=ciHigh; % 按实际收敛样本计算区间并显式报告有效样本数。
results.runtime=runtime; results.iterations=iterations; results.finalChange=finalChange; results.converged=converged; % 输出性能和收敛数据。
results.noiseVariance=noiseVariance; results.estimatedAlpha=estimatedAlpha; results.estimatedTau=estimatedTau; % 输出精度诊断。
results.truncationNMSE=truncation; results.modelError=modelError; results.example=example; % 输出模型误差与实例。
results.delayPower=delayPower; results.meanDelayPower=mean(delayPower,1); % 提供逐信道样本及蒙特卡洛平均时延功率分布。
results.exampleDiagnostics=diagnostics; results.status=sprintf('Converged-only NMSE is primary; valid counts per [SNR,algorithm] are stored. All-sample means are nmseLinearAll. Model assumptions remain explicit.'); % 标明未收敛样本未参与主曲线均值。
fprintf('各SNR/算法收敛样本数（共%d次）：\n',cfg.nMC); disp(array2table(validCount,'VariableNames',matlab.lang.makeValidName(labels),'RowNames',cellstr(num2str(cfg.snrDb(:))))); % 让收敛覆盖率在MATLAB命令行可见。
if cfg.makePlots % 按需生成NMSE-SNR图和真实等效信道热图。
    curveFig=figure('Color','w','Name','OTFS NMSE vs SNR'); plot(cfg.snrDb,results.nmseDb,'-o','LineWidth',1.4); grid on; % 绘制所有算法曲线。
    xlabel('SNR (dB)'); ylabel('NMSE (dB)'); legend(labels,'Location','best','Interpreter','none'); title([cfg.profile,' / ',cfg.pulse,' / ',cfg.snrDefinition],'Interpreter','none'); % 明确标注场景和SNR定义。
    refGrid=example.channel.hDD; refScale=max(abs(refGrid(:))); fullAxis=example.channel.dFullAxis; fullOrder=mod(fullAxis,cfg.N)+1; % 按首个信道计算统一色标并将整周期多普勒列居中。
    channelFig=figure('Color','w','Name','Complete effective DD channel'); imagesc(fullAxis,0:cfg.M-1,20*log10(max(abs(refGrid(:,fullOrder))/max(refScale,realmin),1e-4))); axis xy; colorbar; caxis([-40,0]); % 展示32x16完整离散等效信道。
    xlabel('Doppler tap (centered)'); ylabel('Delay tap'); title('Full 32 x 16 effective DD channel'); % 横向多普勒负坐标表示模N折返后的抽头。
    supportFig=figure('Color','w','Name','Complete DD physical delay support'); imagesc(fullAxis,0:example.channel.lmax,20*log10(max(abs(example.channel.hFullDD(:,fullOrder))/max(refScale,realmin),1e-4))); axis xy; colorbar; caxis([-40,0]); % 放大完整信道的物理非零时延区。
    xlabel('Doppler tap (centered)'); ylabel('Delay tap'); title('All N Doppler bins within the delay support'); % 观察全部远端分数多普勒泄漏。
    windowFig=figure('Color','w','Name','Paper finite-leakage approximation'); imagesc(example.channel.dAxis,0:example.channel.lmax,20*log10(max(abs(example.channel.hWindow)/max(refScale,realmin),1e-4))); axis xy; colorbar; caxis([-40,0]); % 单独显示原论文有限泄漏窗口。
    xlabel('Doppler tap'); ylabel('Delay tap'); title('Truncated leakage used by paperWindow'); % 不再将局部窗口误称作完整信道。
end % 完成可选绘图。
if cfg.saveOutputs % 按需持久保存本次实验输出。
    if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end % 创建结果目录。
    save(fullfile(cfg.outputDir,'nmse_results.mat'),'results','-v7'); % 保存配置、原始样本、均值及诊断，绝不只存图片。
    if cfg.makePlots % 如果创建了图窗则一并保存。
        savefig(curveFig,fullfile(cfg.outputDir,'nmse_vs_snr.fig')); exportgraphics(curveFig,fullfile(cfg.outputDir,'nmse_vs_snr.png'),'Resolution',300); % 保存可编辑曲线和高分辨率图片。
        savefig(channelFig,fullfile(cfg.outputDir,'channel_full_grid.fig')); exportgraphics(channelFig,fullfile(cfg.outputDir,'channel_full_grid.png'),'Resolution',300); % 保存完整32x16 DD网格。
        savefig(supportFig,fullfile(cfg.outputDir,'channel_full_delay_support.fig')); exportgraphics(supportFig,fullfile(cfg.outputDir,'channel_full_delay_support.png'),'Resolution',300); % 保存包含全部16个多普勒抽头的局部放大图。
        savefig(windowFig,fullfile(cfg.outputDir,'channel_paper_window.fig')); exportgraphics(windowFig,fullfile(cfg.outputDir,'channel_paper_window.png'),'Resolution',300); % 保存原文有限泄漏对照图。
    end % 完成图形保存。
end % 完成结果保存。
disp(array2table(results.nmseDb,'VariableNames',matlab.lang.makeValidName(labels),'RowNames',cellstr(num2str(cfg.snrDb(:))))); % 在命令行显示各SNR的NMSE。
end % 结束实验入口，以下均为本文件的局部function。

function summary=compare_six_scenarios(caseMode,overrides) % 在相同表I参数下逐个调用A/B/D与两种波形，并汇总六情景。
assert(~isfield(overrides,'profile') && ~isfield(overrides,'pulse'),'比较模式固定遍历A/B/D和两种波形，不需指定profile或pulse'); % 防止传入的单情景配置被静默覆盖。
profiles={'NTN-TDL-A','NTN-TDL-B','NTN-TDL-D'}; pulses={'ideal','rectangular'}; % 列依次为A/B/D，行依次为双正交/矩形。
if isfield(overrides,'outputDir'), rootDir=overrides.outputDir; else, rootDir=fullfile(pwd,'otfs_profile_comparison'); end % 比较结果的独立根目录。
if isfield(overrides,'makePlots'), makePlots=logical(overrides.makePlots); else, makePlots=true; end % 用户可在无图形环境中关闭绘图。
if isfield(overrides,'saveOutputs'), saveOutputs=logical(overrides.saveOutputs); else, saveOutputs=true; end % 用户可选择只在工作区取得结果。
cases=cell(2,3); channelInfo=cell(2,3); casePaths=cell(2,3); % 分别保存六组全量统计、首信道和结果路径。
for iw=1:2 % 先固定波形，再遍历三种标准信道模型。
    for ip=1:3 % 每个情景调用原有单情景实验入口，共用相同的SNR和默认参数。
        local=overrides; local.profile=profiles{ip}; local.pulse=pulses{iw}; local.fig3Algorithms=true; % 强制仅切换PDP/波形并使用Fig.3六种算法。
        local.outputDir=fullfile(rootDir,[profiles{ip},'_',pulses{iw}]); local.makePlots=false; local.saveOutputs=saveOutputs; local.verbose=false; % 各情景独立保存MAT文件，总览集中绘图。
        fprintf('比较情景 %d/6：%s，%s。\n',(iw-1)*3+ip,profiles{ip},pulses{iw}); % 给出六次仿真的总体进度。
        cases{iw,ip}=run_otfs_nmse(caseMode,local); % 每次调用都会重置同一个随机种子，同一PDP的两种波形共享物理路径。
        channelInfo{iw,ip}=cases{iw,ip}.example.channel; casePaths{iw,ip}=local.outputDir; % 暴露原始抽头、完整信道和每个实验的文件路径。
    end % 完成当前波形下的A/B/D。
end % 完成六情景模拟。
for ip=1:3 % 验证对比两种脉冲时没有意外换掉同一PDP的物理信道。
    chA=channelInfo{1,ip}; chB=channelInfo{2,ip}; % 读取相同随机种子的理想与矩形首实例。
    assert(isequal(chA.gain,chB.gain) && isequal(chA.nuHz,chB.nuHz) && isequal(chA.delayTap,chB.delayTap),'同一PDP的两种波形没有使用相同物理路径'); % 保证波形配对比较。
end % 完成路径一致性检查。
summary=struct('profiles',{profiles},'pulses',{pulses},'cases',{cases},'channelInfo',{channelInfo},'casePaths',{casePaths}); % 六情景原始MATLAB结果均可用summary.cases{行,列}访问。
summary.snrDb=cases{1,1}.cfg.snrDb; summary.labels=cases{1,1}.labels; summary.rootDir=rootDir; summary.caseMode=caseMode; % 保存公共坐标和文件目录。
summary.nmseDb=nan(numel(summary.snrDb),numel(summary.labels),2,3); summary.validCount=zeros(numel(summary.snrDb),numel(summary.labels),2,3); summary.convergedFraction=zeros(numel(summary.snrDb),numel(summary.labels),2,3); % 统一存为[SNR,算法,波形,信道]四维数组，含有效样本与收敛比例。
summary.pilotOverhead=zeros(2,3); summary.quantizedDelayTap=cell(2,3); summary.meanDelayPower=cell(2,3); % 保存导频保护区占比和每个信道的离散数据。
for iw=1:2 % 将六个独立实验的量化指标整理到同一个变量。
    for ip=1:3 % 每列对应一种PDP。
        r=cases{iw,ip}; ch=channelInfo{iw,ip}; % 读取单个情景结果与首信道。
        summary.nmseDb(:,:,iw,ip)=r.nmseDb; summary.validCount(:,:,iw,ip)=r.validCount; summary.convergedFraction(:,:,iw,ip)=r.convergedFraction; summary.quantizedDelayTap{iw,ip}=ch.profileQuantizedDelayTap; summary.meanDelayPower{iw,ip}=r.meanDelayPower; % 不混合六组NMSE的真值与维度，同时暴露每点有效样本数。
        summary.pilotOverhead(iw,ip)=(2*ch.lmax+1)*(4*ch.knu+1)/(r.cfg.M*r.cfg.N); % 按原文单导频保护区公式记录各PDP开销。
    end % 完成当前波形下各模型统计。
end % 完成六情景汇总。
summary.note='Figure-legend algorithms: OMP, BSBL, VSBL, APCSBL, NCP-VSBL, NNCP-VSBL. BSBL uses block-correlated Gaussian EM with adaptive block reorganization; APCSBL uses data-adaptive adjacent pattern coupling. Reference-paper updates are not fully specified in the target paper, so these two are structural EM implementations, not author-code bit-exact replicas. Paired waveform cases share physical paths. Rectangular waveform follows the target paper DD-domain phase model, not a full time-domain receiver.'; % 记录算法来源边界和矩形波形模型范围。
if makePlots % 绘制六个子图以直接查看不同PDP与脉冲对应的原始趋势。
    nmseFigure=figure('Color','w','Name','Six OTFS NMSE scenarios','Position',[100,100,1600,950]); nmseLayout=tiledlayout(nmseFigure,2,3,'TileSpacing','compact','Padding','compact'); % 建立两行三列宽幅曲线图并预留公共图例。
    fullFigure=figure('Color','w','Name','Six complete DD channels'); % 建立同一布局的完整泄漏信道热图。
    for iw=1:2 % 第一行双正交，第二行矩形。
        for ip=1:3 % 三列依次A、B、D。
            r=cases{iw,ip}; ch=channelInfo{iw,ip}; idx=(iw-1)*3+ip; % 从结果中取得当前情景的坐标和信道。
            axesNMSE=nexttile(nmseLayout,idx); plot(axesNMSE,r.cfg.snrDb,r.nmseDb,'-o','LineWidth',1); grid(axesNMSE,'on'); % 每个子图显示同一情景下的全部对比算法。
            xlabel(axesNMSE,'SNR (dB)'); ylabel(axesNMSE,'NMSE (dB)'); title(axesNMSE,[profiles{ip},' / ',pulses{iw}],'Interpreter','none'); % 图题显式标记信道与波形。
            if idx==6, lgd=legend(axesNMSE,summary.labels,'Location','southoutside','Orientation','horizontal','Interpreter','none'); lgd.Layout.Tile='south'; end % 公共图例占用独立布局行，不挤压最后一个子图。
            axisD=ch.dFullAxis; order=mod(axisD,r.cfg.N)+1; gainScale=max(abs(ch.hDD(:))); gainScale=max(gainScale,realmin); % 将16列完整多普勒重新排列为居中坐标并取每幅图的标量幅度基准。
            figure(fullFigure); subplot(2,3,idx); imagesc(axisD,0:r.cfg.M-1,20*log10(max(abs(ch.hDD(:,order))/gainScale,1e-4))); axis xy; caxis([-40,0]); colorbar; % 固定色条以比较六种场景的离散结构。
            xlabel('Doppler tap'); ylabel('Delay tap'); title([profiles{ip},' / ',pulses{iw}],'Interpreter','none'); % 完整网格中的空白时延行也如实保留。
        end % 完成一行三个子图。
    end % 完成两种波形的六个子图。
    pdpFigure=figure('Color','w','Name','Mean delay power of A B D'); % 新建各信道经过时延离散化后的平均边际能量图。
    for ip=1:3 % 同一PDP两个波形在期望功率上应接近，优先用理想波形曲线看时延结构。
        power=cases{1,ip}.meanDelayPower; powerDb=10*log10(max(power/max(max(power),realmin),1e-5)); % 归一化到当前模型最强时延抽头。
        subplot(1,3,ip); stem(0:numel(powerDb)-1,powerDb,'filled'); grid on; ylim([-50,3]); % 保留离散抽头间隔，突出时延空行和远端弱径。
        xlabel('Delay tap'); ylabel('Mean power (dB)'); title(profiles{ip},'Interpreter','none'); % 指明此图是多次信道实现的均值。
    end % 完成三个PDP边际能量图。
end % 完成可选的综合绘图。
if saveOutputs % 以独立文件保存总览与可重复六情景的原始结果。
    if ~exist(rootDir,'dir'), mkdir(rootDir); end % 创建总览目录。
    save(fullfile(rootDir,'six_scenarios.mat'),'summary','-v7'); % MATLAB可直接读取全部六组参数、信道、样本和NMSE。
    if makePlots % 用户开启作图时才输出图片及可编辑的MATLAB图窗。
        savefig(nmseFigure,fullfile(rootDir,'six_nmse.fig')); exportgraphics(nmseFigure,fullfile(rootDir,'six_nmse.png'),'Resolution',300); % 保存六情景NMSE。
        savefig(fullFigure,fullfile(rootDir,'six_channels.fig')); exportgraphics(fullFigure,fullfile(rootDir,'six_channels.png'),'Resolution',300); % 保存完整DD热图。
        savefig(pdpFigure,fullfile(rootDir,'three_delay_pdp.fig')); exportgraphics(pdpFigure,fullfile(rootDir,'three_delay_pdp.png'),'Resolution',300); % 保存三种模型时延PDP。
    end % 完成可编辑图和PNG保存。
end % 完成全部总览文件保存。
end % 结束六情景比较子函数。

function yDD=apply_channel_dd(xDD,ch,cfg) % 对应式(5)的有限求和或完整周期Dirichlet求和；接收延迟为行、多普勒为列。
yDD=complex(zeros(cfg.M,cfg.N)); % 初始化无噪声接收网格。
if strcmpi(cfg.estimationExtent,'complete'), qAxis=0:cfg.N-1; else, qAxis=-cfg.Nhat:cfg.Nhat; end % 完整模式每个路径恰取N个不同的循环多普勒移位。
for i=1:ch.numberPaths % 逐条路径累加，不依赖等效信道估计向量。
    if strcmpi(cfg.pulse,'ideal') % 双正交波形相位为常数。
        phase=exp(-1i*2*pi*ch.delayTap(i)*ch.nuBin(i)/(cfg.M*cfg.N)); % 论文式(6)第一行。
    else % 矩形波形仅保证CE区l>=lmax的式(6)适用，不能用本分支声称完成整帧BER仿真。
        phase=exp(1i*2*pi*((0:cfg.M-1)'-ch.delayTap(i))*ch.nuBin(i)/(cfg.M*cfg.N)); % 论文式(6)第二行在CE区的相位。
    end % 完成波形相位计算。
    for q=qAxis % 原文模式截断泄漏，完整模式用全部N个频率抽头。
        D=mean(exp(1i*2*pi*(0:cfg.N-1)*(q+ch.kappa(i))/cfg.N)); % 有限几何和与式(5)比值严格等价。
        shifted=circshift(xDD,[ch.delayTap(i),ch.k(i)-q]); % 实现x([k-ki+q]N,[l-li]M)。
        yDD=yDD+ch.gain(i)*D*(phase.*shifted); % 累加路径增益、泄漏、相位及移位发送网格。
    end % 完成当前路径的泄漏求和。
end % 完成全部路径传输。
end % 结束DD路径信道函数。

function A=build_measurement(xDD,ch,cfg,estAxis) % 按式(10)构建测量矩阵，并接受完整或有限的显式多普勒轴。
nr=ch.lmax+1; nd=numel(estAxis); A=complex(zeros(nr*nd)); % 导频观测数等于本次估计向量长度。
for jd=1:nd % 遍历观测多普勒索引。
    for jt=1:nr % 遍历观测延迟索引。
        row=jt+(jd-1)*nr; l=cfg.lp+jt-1; k=cfg.kp+estAxis(jd); % 统一使用选定的0基物理多普勒坐标。
        for hd=1:nd % 遍历未知系数的多普勒坐标。
            for ht=1:nr % 遍历未知系数的延迟坐标。
                col=ht+(hd-1)*nr; A(row,col)=xDD(mod(l-(ht-1),cfg.M)+1,mod(k-estAxis(hd),cfg.N)+1); % 论文psi定义，行列均按同一多普勒轴堆叠。
            end % 完成一个多普勒列内的所有时延。
        end % 完成当前观测的整行字典。
    end % 完成当前观测多普勒下的所有延迟。
end % 完成字典矩阵构建。
end % 结束字典构建函数。

function neighbors=make_neighbors(nr,nd,kind) % 原文邻接与物理多普勒邻接两种模式明确区分。
L=nr*nd; center=(1:L)'; % 一维延迟快变索引。
if strcmpi(kind,'paperCircular') % 默认严格保持目标文式(25)循环邻接，不悄悄修改成新方法。
    neighbors=[mod(center-2,L)+1,center,mod(center,L)+1]; % g=1/2/3依次取前一个、自身、后一个一维精度。
elseif strcmpi(kind,'dopplerOpen') % 可选消融：同一延迟上的多普勒开边界邻接。
    t=mod(center-1,nr); d=floor((center-1)/nr); neighbors=[t+max(d-1,0)*nr+1,center,t+min(d+1,nd-1)*nr+1]; % 边界重复自身，相邻延迟不互连。
else % 拒绝未知先验拓扑。
    error('neighborhood应为paperCircular或dopplerOpen'); % 提供明确拓扑选项。
end % 完成邻接选择。
end % 结束邻接构建函数。

function [h,info]=ls_est(y,A,opts) % LS自检基线，不需要信道先验。
h=A\y; info=basic_info(1,true,NaN,NaN); info.finalChange=0; % 使用反斜杠求解，不显式计算逆矩阵。
end % 结束LS算法。

function [h,info]=omp_est(y,A,opts) % 固定原子预算OMP，不访问真实物理路径数或噪声。
L=size(A,2); h=complex(zeros(L,1)); r=y; selected=[]; scale=sqrt(sum(abs(A).^2,1)); % 初始化残差与列范数。
budget=min(opts.ompAtoms,rank(A)); % 原子预算由配置给出；避免超过矩阵秩。
for it=1:budget % 按预算执行贪婪迭代。
    score=abs(A'*r)./max(scale(:),opts.floor); score(selected)=-Inf; % 归一化相关性并屏蔽已经选择的原子。
    [~,j]=max(score); selected(end+1)=j; hs=A(:,selected)\y; r=y-A(:,selected)*hs; % 扩展支持并重新联合LS拟合。
    if norm(r)<=opts.floor*max(norm(y),1), break; end % 数值零残差时提前终止。
end % 完成OMP支持搜索。
if ~isempty(selected), h(selected)=hs; end % 在所选支持外填零。
info=basic_info(numel(selected),true,NaN,NaN); info.finalChange=0; % OMP以固定支持预算结束，不估计噪声精度。
end % 结束OMP算法。

function [h,info]=bsbl_br_legacy_contiguous(y,A,opts) % 旧版错误连续分块，仅保留作审计对照，不参与Fig.3。
L=size(A,2); m=size(A,1); d0=max(1,min(opts.blockLength,L)); edges=1:d0:(L+1); edges(end+1)=L+1; edges=unique(edges); % 初始块对应每个时延tap的多普勒泄漏簇。
edges=sort(edges); rho=0.2; gamma=ones(numel(edges)-1,1); alpha=opts.alpha0; done=false; % 初始化块功率、块内AR相关及噪声精度。
for it=1:opts.maxIter % BSBL的EM更新及周期性块重组。
    C=complex(zeros(L)); % 先验块协方差。
    for b=1:numel(gamma), ix=edges(b):edges(b+1)-1; B=rho.^abs((1:numel(ix))'-(1:numel(ix))); C(ix,ix)=gamma(b)*B; end % 每块采用学习到相关系数的Toeplitz协方差。
    C=C+opts.floor*eye(L); % 防止块先验协方差奇异。
    P=alpha*(A'*A)+C\eye(L); P=(P+P')/2; [R,flag]=chol(P); assert(flag==0,'BSBL后验精度非正定'); S=R\(R'\eye(L)); h=R\(R'\(alpha*A'*y)); S=(S+S')/2; % 高斯E步。
    moment=abs(h).^2+max(real(diag(S)),0); newGamma=gamma; crossSum=0; diagSum=0; % 系数二阶矩和块相关统计量。
    for b=1:numel(gamma), ix=edges(b):edges(b+1)-1; B=rho.^abs((1:numel(ix))'-(1:numel(ix))); Q=S(ix,ix)+h(ix)*h(ix)'; newGamma(b)=max(real(trace(B\Q))/numel(ix),opts.floor); if numel(ix)>1, crossSum=crossSum+sum(abs(diag(Q,1))); diagSum=diagSum+sum(real(diag(Q))); end; end % 更新各非稀疏块尺度。
    newRho=min(0.95,max(0,crossSum/max(diagSum,opts.floor))); newAlpha=m/max(expected_residual(y,A,h,S),opts.floor); change=norm(newGamma-gamma)/max(norm(gamma),opts.floor); gamma=newGamma; rho=newRho; alpha=newAlpha; % 更新块参数和噪声。
    if mod(it,10)==0 && numel(gamma)>1 % 每十轮按相邻块强度重组支持块大小。
        blockEnergy=zeros(numel(gamma),1); for b=1:numel(gamma), blockEnergy(b)=gamma(b); end % 当前块活动度。
        merge=find(abs(diff(log10(max(blockEnergy,opts.floor))))<0.25,1); if ~isempty(merge), edges(merge+1)=[]; gamma(merge)=mean(gamma(merge:merge+1)); gamma(merge+1)=[]; else, split=find(diff(edges)>2,1); if ~isempty(split), mid=floor((edges(split)+edges(split+1))/2); edges=[edges(1:split),mid,edges(split+1:end)]; gamma=[gamma(1:split);gamma(split);gamma(split+1:end)]; end; end % 相近活动块合并，较长块拆分以适应实际簇宽。
    end % 完成一次块重组。
    if it>=opts.minIter && change<opts.tol, done=true; break; end % 收敛检查。
end % 结束BSBL-BR EM。
C=complex(zeros(L)); for b=1:numel(gamma), ix=edges(b):edges(b+1)-1; C(ix,ix)=gamma(b)*rho.^abs((1:numel(ix))'-(1:numel(ix))); end; P=alpha*(A'*A)+C\eye(L); P=(P+P')/2; [R,flag]=chol(P); assert(flag==0,'BSBL最终后验精度非正定'); h=R\(R'\(alpha*A'*y)); info=basic_info(it,done,alpha,NaN); info.blockEdges=edges; info.blockCorrelation=rho; % 以最终模型重新计算后验均值。
end % 结束BSBL-BR入口。

function [h,info]=bsbl_br(y,A,opts) % 分块SBL：同一时延行跨相邻多普勒列组成块，索引匹配h(:)布局。
L=size(A,2); m=size(A,1); nr=opts.nDelay; nd=opts.nDoppler; assert(nr*nd==L,'BSBL的时延/多普勒维度与h不匹配'); % h按delay-fastest顺序堆叠。
edges=[1,nd+1]; gamma=ones(nr,1); rho=0.2; alpha=opts.alpha0; done=false; change=Inf; % 每条时延初始为一个完整多普勒块。
for it=1:opts.maxIter % BSBL EM，按多普勒方向更新相关块。
    C=complex(zeros(L)); % 建立先验协方差。
    for b=1:numel(edges)-1, dop=edges(b):edges(b+1)-1; B=rho.^abs((1:numel(dop))'-(1:numel(dop))); for ell=1:nr, ix=ell+(dop-1)*nr; C(ix,ix)=gamma(ell,b)*B; end; end % 按delay分组，再沿Doppler建立Toeplitz块。
    C=C+opts.floor*eye(L); P=alpha*(A'*A)+C\eye(L); P=(P+P')/2; [R,flag]=chol(P); assert(flag==0,'BSBL后验精度非正定'); S=R\(R'\eye(L)); h=R\(R'\(alpha*A'*y)); S=(S+S')/2; % E步。
    newGamma=gamma; crossSum=0; diagSum=0; % 块二阶矩与块内相关统计量。
    for b=1:size(gamma,2), dop=edges(b):edges(b+1)-1; B=rho.^abs((1:numel(dop))'-(1:numel(dop))); for ell=1:nr, ix=ell+(dop-1)*nr; Q=S(ix,ix)+h(ix)*h(ix)'; newGamma(ell,b)=max(real(trace(B\Q))/numel(dop),opts.floor); if numel(dop)>1, crossSum=crossSum+sum(abs(diag(Q,1))); diagSum=diagSum+sum(real(diag(Q))); end; end; end % 每个时延行、多普勒段单独估计块功率。
    newRho=min(0.95,max(0,crossSum/max(diagSum,opts.floor))); newAlpha=m/max(expected_residual(y,A,h,S),opts.floor); change=max([norm(newGamma-gamma)/max(norm(gamma),opts.floor),abs(newRho-rho),abs(newAlpha-alpha)/max(alpha,opts.floor)]); gamma=newGamma; rho=newRho; alpha=newAlpha; % M步，并把相关系数和噪声精度纳入停止准则。
    if mod(it,10)==0 % 根据沿多普勒方向的平均能量重组块。
        segEnergy=mean(gamma,1); merge=find(abs(diff(log10(max(segEnergy,opts.floor))))<0.25,1); if ~isempty(merge), edges(merge+1)=[]; gamma(:,merge)=mean(gamma(:,merge:merge+1),2); gamma(:,merge+1)=[]; elseif rho<0.35, widths=diff(edges); [wid,b]=max(widths); if wid>2, mid=floor((edges(b)+edges(b+1))/2); edges=[edges(1:b),mid,edges(b+1:end)]; gamma=[gamma(:,1:b),gamma(:,b),gamma(:,b+1:end)]; end; end % 块边界定义在Doppler轴上，不混合相邻时延。
    end % 完成本轮可能的块重组。
    if it>=opts.minIter && change<opts.tol, done=true; break; end % 检查关键超参数变化。
end % 结束BSBL迭代。
C=complex(zeros(L)); for b=1:size(gamma,2), dop=edges(b):edges(b+1)-1; B=rho.^abs((1:numel(dop))'-(1:numel(dop))); for ell=1:nr, ix=ell+(dop-1)*nr; C(ix,ix)=gamma(ell,b)*B; end; end; C=C+opts.floor*eye(L); P=alpha*(A'*A)+C\eye(L); P=(P+P')/2; [R,flag]=chol(P); assert(flag==0,'BSBL最终后验精度非正定'); h=R\(R'\(alpha*A'*y)); info=basic_info(it,done,alpha,NaN); info.finalChange=change; info.blockEdges=edges; info.blockCorrelation=rho; % 保存最终收敛指标和多普勒块边界。
end % 结束BSBL新实现。

function [h,info]=apcsbl_adaptive(y,A,opts) % 自适应邻接耦合SBL，EM更新局部精度及数据驱动的双侧耦合。
L=size(A,2); m=size(A,1); alphaLocal=ones(L,1); betaLeft=0.5*ones(L,1); betaRight=0.5*ones(L,1); alpha=opts.alpha0; done=false; omega=1./max(alphaLocal,opts.floor); [h,S]=gaussian_posterior(y,A,alpha,omega); % 初始化耦合精度与高斯后验。
for it=1:opts.maxIter % APCSBL的交替后验和自适应超参数更新。
    moment=abs(h).^2+max(real(diag(S)),0); left=[moment(1);moment(1:end-1)]; right=[moment(2:end);moment(end)]; % 邻接后验能量。
    crossL=abs([0;conj(h(1:end-1)).*h(2:end)]); crossR=abs([conj(h(1:end-1)).*h(2:end);0]); % 邻接系数相关幅度。
    scaleL=sqrt(max(moment.*left,opts.floor)); scaleR=sqrt(max(moment.*right,opts.floor)); newLeft=min(0.9,max(0,crossL./scaleL)); newRight=min(0.9,max(0,crossR./scaleR)); % 按后验相邻相关性自适应耦合强度。
    neighborL=[alphaLocal(1);alphaLocal(1:end-1)]; neighborR=[alphaLocal(2:end);alphaLocal(end)]; coupledEnergy=moment+newLeft.*left+newRight.*right; newLocal=1./max(coupledEnergy,opts.floor); newOmega=newLocal+newLeft.*neighborL+newRight.*neighborR; % pattern-coupled Gaussian先验精度。
    newAlpha=m/max(expected_residual(y,A,h,S),opts.floor); change=norm(newOmega-omega)/max(norm(omega),opts.floor); alphaLocal=newLocal; betaLeft=newLeft; betaRight=newRight; omega=max(newOmega,opts.floor); alpha=newAlpha; [h,S]=gaussian_posterior(y,A,alpha,omega); % EM更新并重新估计后验。
    if it>=opts.minIter && change<opts.tol, done=true; break; end % 收敛判定。
end % 完成自适应模式耦合更新。
info=basic_info(it,done,alpha,NaN); info.finalChange=change; info.betaLeft=betaLeft; info.betaRight=betaRight; info.omega=omega; % 保存数据自适应耦合参数及停止误差。
end % 结束APCSBL入口。

function [h,info]=sbl_em(y,A,opts) % 独立高斯方差的标准复数SBL-EM基线。
L=size(A,2); m=size(A,1); variance=ones(L,1)/opts.gamma0; alpha=opts.alpha0; done=false; % 初始化局部方差和噪声精度。
for it=1:opts.maxIter % 交替E步与M步。
    [h,S]=gaussian_posterior(y,A,alpha,1./variance); moment=abs(h).^2+real(diag(S)); % E步：均值、全协方差及二阶矩。
    newVariance=max(moment,opts.floor); newAlpha=m/max(expected_residual(y,A,h,S),opts.floor); % M步：局部方差及复高斯噪声精度更新。
    change=norm(newVariance-variance)/max(norm(variance),opts.floor); variance=newVariance; alpha=newAlpha; % 记录相对改变并更新参数。
    if it>=opts.minIter && change<opts.tol, done=true; break; end % 统一迭代阈值，保留收敛标记。
end % 完成SBL-EM循环。
[h,~]=gaussian_posterior(y,A,alpha,1./variance); info=basic_info(it,done,alpha,NaN); info.finalChange=change; % 使用最终参数重新计算输出均值并保存停止误差。
end % 结束独立SBL-EM。

function [h,info]=vb_sbl(y,A,opts) % 独立Gamma-Gaussian复数VB基线，不等同Wang文中的联合数据检测器。
L=size(A,2); m=size(A,1); a=1e-10; b=a; c=a; d=a; gamma=opts.gamma0*ones(L,1); alpha=opts.alpha0; done=false; % 采用目标论文脚注同量级超先验。
[h,S]=gaussian_posterior(y,A,alpha,gamma); % 初始化先验和似然匹配的高斯后验。
for it=1:opts.maxIter % 标准独立均值场更新。
    alpha=(c+m)/max(d+expected_residual(y,A,h,S),opts.floor); % 噪声Gamma形状率更新。
    [h,S]=gaussian_posterior(y,A,alpha,gamma); moment=abs(h).^2+real(diag(S)); % 高斯信道后验更新。
    newGamma=(a+1)./max(b+moment,opts.floor); change=norm(newGamma-gamma)/max(norm(gamma),opts.floor); gamma=newGamma; % 每个复高斯贡献一个形状计数。
    if it>=opts.minIter && change<opts.tol, done=true; break; end % 满足统一阈值则停止。
end % 完成独立VB循环。
[h,~]=gaussian_posterior(y,A,alpha,gamma); info=basic_info(it,done,alpha,NaN); info.finalChange=change; % 输出最终参数对应的后验均值和停止误差。
end % 结束独立VB算法。

function [h,info]=ncp_vb(y,A,opts) % Dai2019作者公开先验的固定网格单测量向量版本。
[h,info]=pattern_vb_core(y,A,opts,false); % 关闭全局tau和alpha-tau耦合；不包含原作者角度网格细化。
end % 结束NCP固定网格算法。

function [h,info]=nncp_vsbl(y,A,opts) % 目标论文NNCP-VSBL：严格采用式(21)-(36)的更新顺序。
[h,info]=pattern_vb_core(y,A,opts,true); % 开启p(tau|alpha)=Ga(e,alpha)和全局精度。
end % 结束目标论文算法入口。

function [h,info]=pattern_vb_core(y,A,opts,useCoupling) % 共享数学运算，两个独立入口明确区分两种联合分布。
L=size(A,2); m=size(A,1); nb=opts.neighbors; a=1e-10; b=a; c=a; d=a; e=a; % 超先验采用论文脚注参数。
gamma=opts.gamma0*ones(L,1); phi=ones(L,3)/3; alpha=opts.alpha0; tau=1; done=false; % 式(13)的三选一先验均匀初始化。
if useCoupling, tau=opts.tau0; end % 仅NNCP具有全局精度变量。
omega=tau*sum(phi.*gamma(nb),2); [h,S]=gaussian_posterior(y,A,alpha,omega); % 初始Omega和后验，避免使用真实噪声初始化。
traceAlpha=zeros(opts.maxIter,1); traceTau=traceAlpha; traceChange=traceAlpha; % 保存检查用迭代轨迹。
for it=1:opts.maxIter % 按目标论文算法1执行参数块更新。
    oldOmega=omega; oldTau=tau; % 式(28)需要旧tau，停止准则需要旧Omega。
    if useCoupling % NNCP的条件Gamma先验给噪声后验额外计数e和率tau。
        alpha=(c+e+m)/max(d+tau+expected_residual(y,A,h,S),opts.floor); % 原文式(21)-(22)，m=L时完全对应原式。
    else % NCP没有全局tau的条件先验。
        alpha=(c+m)/max(d+expected_residual(y,A,h,S),opts.floor); % Dai公开代码的单测量向量噪声更新。
    end % 完成噪声精度更新。
    [h,S]=gaussian_posterior(y,A,alpha,omega); moment=abs(h).^2+real(diag(S)); % 原文式(23)-(25)，moment对应varpi。
    shape=a*ones(L,1); rate=b*ones(L,1); % Gamma局部精度的形状和率。
    for g=1:3 % 按“系数j选择精度nb(j,g)”反向收集贡献，避免左右移位写反。
        shape=shape+accumarray(nb(:,g),phi(:,g),[L,1]); % 原文式(27)，重复边界也能正确累计。
        rate=rate+oldTau*accumarray(nb(:,g),phi(:,g).*moment,[L,1]); % 原文式(28)，使用旧tau和本轮varpi。
    end % 完成三个选择对局部精度的贡献。
    gamma=shape./max(rate,opts.floor); expectedLogGamma=psi(shape)-log(max(rate,opts.floor)); % 原文式(29)-(30)，不能把Eloggamma替换为logEgamma。
    if useCoupling % 更新全局Gamma精度。
        shapeTau=e+sum(phi(:)); rateTau=alpha+sum(sum(phi.*gamma(nb).*moment)); tau=shapeTau/max(rateTau,opts.floor); % 原文式(31)-(34)，sum(phi(:))=L。
        expectedLogTau=psi(shapeTau)-log(max(rateTau,opts.floor)); % 原文式(35)，在三个模式的softmax中是共同项。
    else % NCP等价固定tau=1而不推断该变量。
        tau=1; expectedLogTau=0; % 不错误保留alpha-tau耦合。
    end % 完成全局精度更新。
    score=expectedLogGamma(nb)+expectedLogTau-tau*gamma(nb).*moment; % 原文式(36)的对数消息。
    score=score-max(score,[],2); phi=exp(score); phi=phi./sum(phi,2); % 行softmax防止指数上溢；均匀先验常数抵消。
    omega=tau*sum(phi.*gamma(nb),2); change=max(abs(omega-oldOmega))/max(max(abs(omega)),opts.floor); % 原文对角Omega的谱范数相对停止准则。
    traceAlpha(it)=alpha; traceTau(it)=tau; traceChange(it)=change; % 保存每一轮精度与收敛误差。
    if it>=opts.minIter && change<opts.tol, done=true; break; end % 原文阈值及安全最小迭代数。
end % 完成NNCP/NCP变分循环。
info=basic_info(it,done,alpha,tau); info.finalChange=change; info.phi=phi; info.gamma=gamma; info.omega=omega; % 保存先验状态和最终停止误差以便人工核对。
info.traceAlpha=traceAlpha(1:it); info.traceTau=traceTau(1:it); info.traceChange=traceChange(1:it); % 删除预分配数组中尚未运行的尾部。
% 保持论文算法1输出当前q(h)均值；不额外添加一次未在迭代顺序中的后验更新。
end % 结束邻域耦合VB核心。

function [mu,S]=gaussian_posterior(y,A,alpha,omega) % 求式(24)-(25)的复高斯后验，避免inv。
L=size(A,2); gram=A'*A; rhs=alpha*(A'*y); offdiag=gram-diag(diag(gram)); % 计算Hermitian Gram矩阵和右端项。
if norm(offdiag,'fro')<=1e-13*max(norm(gram,'fro'),1) % 单导频充分保护时Gram严格对角，这是数学等价加速。
    diagS=1./(alpha*real(diag(gram))+omega); S=diag(diagS); mu=diagS.*rhs; % 直接按分量求解，与矩阵求逆公式完全一致。
else % 一般字典使用Cholesky解线性方程。
    precision=alpha*gram+diag(omega); precision=(precision+precision')/2; [R,flag]=chol(precision); % 只消除舍入导致的非Hermitian误差。
    assert(flag==0,'后验精度矩阵非正定，检查精度或输入矩阵'); S=R\(R'\eye(L)); mu=R\(R'\rhs); % 不静默添加正则化改变原算法。
end % 完成分支求解。
S=(S+S')/2; % 保持数值Hermitian对称。
end % 结束高斯后验求解函数。

function energy=expected_residual(y,A,h,S) % 噪声Gamma后验需要的期望残差能量。
energy=real(norm(y-A*h)^2+sum(sum(conj(A).*(A*S)))); % 等于||y-Ah||^2+tr(A*S*A^H)，不能漏掉协方差项。
assert(isfinite(energy) && energy>=-1e-9,'期望残差非法'); energy=max(energy,0); % 仅去除数值舍入产生的极小负数。
end % 结束期望残差计算。

function info=basic_info(iterations,converged,alpha,tau) % 统一所有算法的诊断接口。
info=struct('iterations',iterations,'converged',converged,'alpha',alpha,'tau',tau); % 无此参数的基线以NaN标记而非伪造估计值。
end % 结束诊断结构构建。
