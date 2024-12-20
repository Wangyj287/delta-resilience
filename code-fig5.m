

%% subfigure a
clc, clear, close all
color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529];
color_mean = color_all(4,:);
color_between = color_all(5,:);

jj=0;
z1=xlsread('../西北江站点_YH.xlsx','马口');
z_mk=z1(133:end,3); % MK
for i=1:12:length(z_mk)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_mk(jj)=mean(z_mk(ii1));
end
jj=0;
b1=xlsread('../西北江站点_YH.xlsx','三灶');
z_sz=b1(133:end,3); % SZ
for i=1:12:length(z_sz)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_sz(jj)=mean(z_sz(ii1));
end
S1=(m_mk-m_sz)/138500; % MK-SZ
S1=S1';

jj=0;
z2=xlsread('../西北江站点_YH.xlsx','三水');
Z_ss=z2(133:end,3); % SS
for i=1:12:length(Z_ss)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_ss(jj)=mean(Z_ss(ii1));
end
jj=0;
z3=xlsread('../西北江站点_YH.xlsx','黄埔');
Z_fbc=z3(133:end,3);
for i=1:12:length(Z_fbc)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_fbc(jj)=mean(Z_fbc(ii1));
end
S4=(m_ss-m_fbc)/61160; % SS-FBC
S4=S4';

jj=0;
z3=xlsread('../西北江站点_YH.xlsx','横门');
Z_hm=z3(133:end,3); % NS
for i=1:12:length(Z_hm)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_hm(jj)=mean(Z_hm(ii1));
end
S2=(m_ss-m_hm)/105000; % SS-HM
S2=S2';

jj=0;
z4=xlsread('../西北江站点_YH.xlsx','蚬沙（南华）');
Z_nh=z4(133:end,3); % NH
for i=1:12:length(Z_nh)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_nh(jj)=mean(Z_nh(ii1));
end
jj=0;
z5=xlsread('../西北江站点_YH.xlsx','南沙');
m3=find(z5(:,1)==1966);
Z_ns=z5(m3:end,3); % NS
for i=1:12:length(Z_ns)
    i1=i; jj=jj+1;   i2=i1+11; ii1=i1:i2;
    m_ns(jj)=mean(Z_ns(ii1));
end
S3=(m_nh-m_ns)/60999; % NH-NS
S3=S3';

MK_w=S1;
MK_w=(MK_w-min(MK_w))./(max(MK_w)-min(MK_w));
SS_w=S2;
SS_w=(SS_w-min(SS_w))./(max(SS_w)-min(SS_w));
SS_v=S3;
SS_v=(SS_v-min(SS_v))./(max(SS_v)-min(SS_v));
SS_u=S4;
SS_u=(SS_u-min(SS_u))./(max(SS_u)-min(SS_u));
dd1=((1989-1966)+0); 
dd2=((1999-1966)+0); 
dd3=((1992-1966)+0); 
dd4=((2010-1966)+0);
cum_S1=cumsum(MK_w-mean(MK_w),1);
cum_S2=cumsum(SS_w-mean(SS_w),1);
cum_S3=cumsum(SS_v-mean(SS_v),1);
cum_S4=cumsum(SS_u-mean(SS_u),1);

a=xlsread('月均流量与水位.xlsx',1);
b=xlsread('月均流量与水位.xlsx',2);
Z=b(73:end,4:2:end);
jj=0;
for i=1966:2016
    for j=1:12
        jj=jj+1;
        time(jj)=datenum(i,j,1);
    end
end
S=(Z(:,2)-Z(:,1))./42367;
s=(S-min(S))./(max(S)-min(S));
S1=mean(reshape(S,12,[]));
s1=(S1-min(S1))./(max(S1)-min(S1));
cum_s=cumsum(s);
ds1=s1-mean(s1);
cum_ds1=cumsum(ds1);

figure;

xlim([0 52]);
set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman','FontSize', 14);
set(gca, 'XTicklabel', {'1966','1970','1980','1990','2000','2010','2016'}, 'Xtick', [1 5 15 25 35 45 51]);
xlabel('\fontname{Times New Roman}Year');
ylabel('\fontname{Times New Roman}The cumulative departure of \newline      water level gradient','FontSize', 16);
hold on;
ylim([-2 6]);
hold on;
combined_array = [cum_S1  cum_S2  cum_S3 cum_S4  cum_ds1'];
average_combined = mean(combined_array, 2);
p1=plot(average_combined, '-.bs');
std_combined = 0.5*std(average_combined); 
upper_combined = average_combined + std_combined;
lower_combined = average_combined - std_combined;
h1=fill([1:length(average_combined) fliplr(1:length(average_combined))], ...
     [upper_combined; flipud(lower_combined)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
set(h1,'Facecolor',color_mean,'FaceAlpha',0.3,'EdgeColor','none');
dd1990 = (1990 - 1966) + 1;  
xline(dd1990, '--k', 'LineWidth', 2);
fill([0 dd1990 dd1990 0], [-2 -2 6 6], color_all(5,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([dd1990 length(average_combined)+1 length(average_combined)+1 dd1990], [-2 -2 6 6], color_all(4,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
text(0.914411696449289, -0.9,'\fontname{Times New Roman}Periods with \newlineless human interventions','FontSize', 14);
text(27, -0.9,'\fontname{Times New Roman}Periods with \newlinesignificant human interventions','FontSize', 14);
yyaxis right;
AR1=[-0.153150780374404,-0.303677304570508,-0.501007799961046,-0.588784143394729,-0.677525242064461,-0.614345190594898,-0.730816384362160,-0.838163826339333,-0.907002935625822,-0.950581847558853,-0.897213129875108,-0.938613937620038,-1.00413980031265,-1.16826840914145,-1.09507426659856,-1.00974700710036,-1.05835540793668,-1.30887124247190,-1.43956433818701,-1.42776591188674,-1.38149692466286,-1.39776120135547,-1.44572975270128,-1.47183034700897,-1.38617591292189,-1.13714627357726,-0.912888253107720,-0.738136718362528,-0.843937975008094,-0.845490066622655,-0.771743210667847,-0.689686038984160,-0.692388258987678,-0.645716892883255,-0.558260529287826,-0.488153789410292,-0.518630523640631,-0.554691897710551,-0.485812683262693,-0.354871756203158,-0.214634792524912,-0.222793299412930,-0.0706635219803027,0.0562825997953821,0.204120163424156,0.131945948438435,0.0621814652269670,3.92543140849609e-16]
windowSize = 6;
AR1 = movmean(AR1, windowSize);
std_AR1 =  0.5*std(AR1);
AR1_extended = [NaN NaN AR1 NaN];
validIndices = ~isnan(AR1_extended);
x_valid = 1:length(AR1_extended);
x_valid = x_valid(validIndices);
AR1_plot = AR1_extended(validIndices);
std_AR1 = std(AR1_plot, 'omitnan');  
upper_AR1 = AR1_plot + std_AR1;
lower_AR1 = AR1_plot - std_AR1;
hold on;
yyaxis right;
p2 = plot(x_valid, AR1_plot, '--or', 'DisplayName', 'Mean AR1');
h2=fill([x_valid, fliplr(x_valid)], ...
     [upper_AR1, fliplr(lower_AR1)], ...
     'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
set(h2,'Facecolor',color_mean,'FaceAlpha',0.3,'EdgeColor','none'); 
xlim([0 52]);
ylabel('\fontname{Times New Roman}The cumulative departure of AR1','FontSize', 14);
legend([p1, p2],{'Mean water level gradient','Mean AR1'});
legend boxoff;
set(gcf,"Position",[10 10 800 400]);
print(gcf, '累积曲线对比', '-dtiff', '-r300')

