function gr_clue(fn_settings,fn_data)
%
% Assessment of the discrete size range from experimental data by the metric proposed by
% Olenius et al., Sci. Rep. 8, 14160, doi:10.1038/s41598-018-32610-z (2018)
% Growth rate (GR) calculations are not valid below the assessed threshold sizes, and these small sizes cannot be described by
% condensational growth/aerosol dynamics models (instead they require molecule-by-molecule growth models such as ACDC)
%
% Input:
%    fn_settings:   name of the .m file that contains variables for input settings
%    fn_data:       name of the .m file that loads (and possibly pre-processes) particle size distribution data
%
% The input .m files are run inside this function, i.e. all variables given by/calculated within them are available here
% The example input files contain more info on the settings and data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n-----------------------------------------------------------------------------------------------')
fprintf('\n                                            GR-CLUE                                            ')
fprintf('\ntool to assess validity of Growth Rate determined from molecular CLUster distribution Evolution')
fprintf('\n-----------------------------------------------------------------------------------------------\n')
fprintf('\n- Below the obtained approximate threshold size(s), growth rate (GR) calculations are not expected to be valid')
fprintf('\n- The sizes depend on ambient conditions such as vapor concentrations, temperature, sink, ...')
fprintf('\n- For more info, run ''help gr_clue''\n\n')
fprintf('\n-----------------------------------------------------------------------------------------------\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default settings

l_show_fits=1;
l_pause=0;

% Size bins considered in the analysis
nbin_min=1;
nbin_max=nan;

% Bins for which the metric is assessed
nbin_start=1;
nbin_end=nan;

% Bin appearance times
apptime_bin=nan;

% Possible fit weight factors
array_weight=nan;
extra_weight=nan;
nadj_weight=nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the input settings

run(fn_settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the input particle size distribution data

run(fn_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add/derive some more parameters

get_constants;

% Assumed number of molecule types (for average properties)
ntypes=length(mv);

% Number of bins and bin indices
ibins=size(dNdlogDp,2);
if isnan(nbin_max), nbin_max=ibins; end
if isnan(nbin_end), nbin_end=ibins; end

% Option for a single size distribution without time dependence (typically a steady state)
if size(dNdlogDp,1) == 1
    time_vector=[0];
    time_tag='last';
end
if size(time_vector,1) > 0
    time_vector=time_vector';
end

% Total concentrations in bins (for e.g. test plots)
%N_conc=bsxfun(@times,dNdlogDp,diff(log10(Dp_limits_bin)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mobility-mass diameter conversions
% The input data is assumed to apply mobility diameters, while mass diameters
% are used for calculations

% Use a simple relation
Dp_mass=Dp_mean_bin-0.3*1e-9;                   % Mass diameter; m
Dp_limits_mass=Dp_limits_bin-0.3*1e-9;

Dp_mass_mon=sum((6/pi/rho*mv).^(1/3))/ntypes;   % Mean molecular mass diameter

% Variables for converting the derivatives of the size distribution
dd0=(Dp_mass_mon^3./(3*Dp_mass.^2))';           % dDp_mass/dn
dd1=(-2*Dp_mass_mon^3./(3*Dp_mass.^3))';        % d(dDp_mass/dn)/dDp_mass
dd2=(2*Dp_mass_mon^3./(Dp_mass.^4))';           % d2(dDp_mass/dn)/dDp_mass2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Criterion for "ok" R2 value (based on how the fits look like)
r2_th=0.90;
l_r2_low_str=['Less good fit (R^2 < ',num2str(r2_th),')'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure settings

Fpos=[200 400 800 500];

Lwidth=1.5;
Msize=10;
Fsize=14;

Mstyle_point='x'; % For marking specific points on lines/colormaps
Msize_point=16;

colors_for_runs=[0 0.45 0.75; 0 0.75 0.75; 1 0.7 0.15; 1 0 0; 1 0 1];

legend_for_bins=strcat(cellfun(@(x) num2str(x,'%.1f'),num2cell(Dp_limits_bin(1:end-1)*1e9),'UniformOutput',false),' - ',...
    cellfun(@(x) num2str(x,'%.1f'),num2cell(Dp_limits_bin(2:end)*1e9),'UniformOutput',false),' nm');

% Lines to indicate different thresholds (e.g. 1%, 5%, ...)
fig_th=sort(unique([1e-2, 5e-2, small_th]),'ascend');
legend_for_th='';
for nth=1:length(fig_th)
    legend_for_th=[legend_for_th, num2str(fig_th(nth))];
    if nth < length(fig_th)
        legend_for_th=[legend_for_th, ', '];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Surface plot (time, diameter, conc.)

if length(time_vector) > 1

    figure(10)
    set(gcf,'Position',Fpos)
    set(gcf,'Color','white')
    box on
    hold on
    pcolor(time_vector,Dp_limits_bin(1:end-1)*1e9,log10(dNdlogDp'))
    %pcolor(time_vector,Dp_limits_bin(1:end-1)*1e9,log10(N_conc'))
    set(gca,'LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize);
    set(gca,'YScale','log')
    set(gca,'Layer','top')
    colormap(jet)
    shading flat

    caxis(lim_for_caxis);
    xlabel('{\itt} (min)')
    xlim([time_vector(1) time_vector(end)])
    ylabel('{\itd}_{p,mob} (nm)')
    ylim([Dp_limits_bin(nbin_min)*1e9 Dp_limits_bin(nbin_max+1)*1e9])
    %set(gca,'YTick',1:10)

    cb=colorbar('LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize);
    ylabel(cb,'d{\itC}/dlog {\itd}_{p,mob} (cm^{-3})','FontWeight','normal','FontSize',Fsize);
    %ylabel(cb,'{\itC} (cm^{-3})','FontWeight','normal','FontSize',Fsize);
    set(cb,'YTick',lim_for_caxis(1):lim_for_caxis(end))
    set(cb,'YTickLabel',cellfun(@(x) num2str(x,'10^{%d}'),num2cell(lim_for_caxis(1):lim_for_caxis(end)),'UniformOutput',false))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bin appearance times, if not given as input: very rough estimate determined by the time resolution
% (i.e. no interpolation over time to increase the time resolution)

if strcmp(time_tag,'tapp')

    %apptime_index_bin=nan(1,ibins);

    if all(isnan(apptime_bin))

        disp('! No input given for bin appearance times - making a rough assessment')
        disp(['  Determining t_app at ',num2str(referenceC_fraction*100),'% of max. conc.'])

        maxC_bin=max(dNdlogDp,[],'omitnan');
        initC_bin=dNdlogDp(1,:);
        referenceC_bin=initC_bin+referenceC_fraction*(maxC_bin-initC_bin);

        % Appearance times for each size bin
        for nbin=1:ibins
            for nt=1:length(time_vector)
                if dNdlogDp(nt,nbin) >= referenceC_bin(nbin)
                    apptime_bin(nbin)=time_vector(nt);
                    %apptime_index_bin(nbin)=nt;
                    break
                end
            end
        end

    % else
    %     
    %     for nbin=1:ibins
    %         if ~isnan(apptime_bin(nbin))
    %             [~,nt]=min(abs(time_vector-apptime_bin(nbin)));
    %             apptime_index_bin(nbin)=nt;
    %         end
    %     end

    end
    
else
    apptime_bin=nan(1,ibins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Particle distribution and its derivatives at given times

if l_show_fits == 0
    fprintf('\nNot showing the fits - remember to ensure that they look decent\n')
end

fit_xaxis=Dp_mass';
% Determine the number density via Dp as DeltaC/DeltaDp (e.g. synthetic molecular data)
%Deltax=diff(Dp_limits_mass)';

nbin_ind=nbin_start:nbin_end;

% Times at which the derivatives are determined
if strcmp(time_tag,'tapp')
    times_for_ddDp_ratio=apptime_bin(nbin_ind);
else
    if strcmp(time_tag,'all')
        times_for_ddDp_ratio=time_vector;
    elseif strcmp(time_tag,'last')
        times_for_ddDp_ratio=time_vector(end);
    end
    nbin_ind_tmp=nbin_ind;
end
nt_tot=length(times_for_ddDp_ratio);

% The metric
ddDp_ratio_Dp2molec=nan(1,ibins);
ddDp_ratio_Dp2molec_vs_time=nan(nt_tot,ibins);

% Some goodness-of-fit statistics (rough)
r2=nan(1,ibins);
r2_vs_time=nan(nt_tot,ibins);

ibins_incl_min=Inf;
ibins_incl_max=0;
ibins_incl_def=10;

for nt=1:nt_tot
    
    % Determine the derivatives at time times_for_ddDp_ratio(nt)
    tval=times_for_ddDp_ratio(nt);
    [tdiffval,tind]=min(abs(time_vector-tval));
    
    % Go through the bins to apply a piecewise fit
    if strcmp(time_tag,'tapp')
        nbin_ind_tmp=nbin_ind(nt);
    end
    
    for nbin=nbin_ind_tmp
        
        %conc_tmp=N_conc(tind,:);
        %conc_tmp=conc_tmp./Deltax'; % Convert to number per size interval; cm^-3 / molec. or m
        conc_tmp=dNdlogDp(tind,:)./(log(10)*Dp_mean_bin); % Convert dNdlogDp to dNdDp

        lfit_ok=0;
        ddDp_tmp=nan(ibins,1);
        d2dDp2_tmp=nan(ibins,1);
        ddDp_ratio_Dp2molec_tmp=nan(ibins,1);
        
        % Fit a function on the distribution
        
        conc_tmp(isnan(conc_tmp))=0;
        
        %%%% 3rd degree polynomial %%%%
        
        % Some ballpark choices for number of bins included in the fit
        % The function is applied around one bin at a time to obtain better fit
        if length(conc_tmp(conc_tmp>0)) >= 100 % Set double amount for very large number of bins (unrealistic, though)
            nadj=ibins_incl_def;
        else
            nadj=ceil(ibins_incl_def/2);
        end
        
        % Use in the fit nadj adjacent points on both sides of nbin, if possible
        fwn1=max(nbin-nadj,nbin_min);
        fwn2=min(fwn1+2*nadj,nbin_max);
        fwn1=max(fwn2-2*nadj,nbin_min);

        % Exclude empty/undefined bins
        zero_points=find(conc_tmp==0);
        excluded_points=setdiff(1:ibins,fwn1:fwn2);
        excluded_points=sort(unique([excluded_points, zero_points]),'ascend');
        
        % Skip the analysed bin if there is not enough data
        ibins_incl=length(conc_tmp)-length(excluded_points);
        
        % Skipped cases:
        % Only up to 4 data points (can be trivially fitted on a 3rd degree polynomial)
        % Current bin and/or adjacent bins missing
        if ibins_incl < 5 || ismember(nbin,excluded_points) || isequal([nbin-1,nbin+1],intersect([nbin-1,nbin+1],excluded_points))
            continue
        end
        ibins_incl_min=min(ibins_incl,ibins_incl_min);
        ibins_incl_max=max(ibins_incl,ibins_incl_max);
        
        % Possible weights
        if isnan(array_weight)
            fweights=ones(1,ibins);
        else
            fweights=array_weight;
        end
        % Give more weight to data points close to the current bin
        if ~isnan(extra_weight)
            ind_tmp=max(nbin-nadj_weight,nbin_min):min(nbin+nadj_weight,nbin_max);
            fweights(ind_tmp)=extra_weight*fweights(ind_tmp);
        end

        [fit_log10,gof]=fit(fit_xaxis,log10(conc_tmp'),'poly3','weights',fweights,'exclude',excluded_points,'Normalize','on');
        lfit_ok=1;

        %%%% Take the derivatives of the fit %%%%
        
        % Analytical derivatives
        [ddDp_log10,d2dDp2_log10]=differentiate(fit_log10,fit_xaxis); % log10(cm^3)/bin or nm, log10(cm^3)/bin or nm^2
        fit_tmp=10.^fit_log10(fit_xaxis);
        ddDp_tmp=log(10)*fit_tmp.*ddDp_log10; % log10(C) -> C
        d2dDp2_tmp=log(10)*fit_tmp.*(d2dDp2_log10+log(10)*ddDp_log10.^2);

        % Convert between Dp and molec. units
        ddDp_ratio_Dp2molec_tmp=dd0.*((dd2.*fit_tmp+2*dd1.*ddDp_tmp+dd0.*d2dDp2_tmp)./(dd1.*fit_tmp+dd0.*ddDp_tmp))+dd1;
        
        if l_show_fits == 1
        
            figure(11)
            clf(11)
            set(gcf,'Position',Fpos)
            set(gca,'LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize)
            set(gcf,'Color','white')
            set(gca,'YScale','log')
            box on
            hold on

            % Plot the distribution
            plot(Dp_mean_bin.*1e9,conc_tmp,'Marker',Mstyle,'MarkerSize',Msize-4,'Color',colors_for_runs(nrun,:),...
                'LineStyle','none','LineWidth',Lwidth,'DisplayName','Data')
            if lfit_ok==1
                Dp_tmp=fit_xaxis(fwn1):1e-12:fit_xaxis(fwn2);
                fit_tmp=10.^fit_log10(Dp_tmp);
                plot((Dp_mean_bin(fwn1):1e-12:Dp_mean_bin(fwn2)).*1e9,fit_tmp,'-k','HandleVisibility','off')

                Lstyle_tmp='ko';
                fit_tmp=10.^fit_log10(fit_xaxis);
                plot(Dp_mean_bin(fwn1:fwn2).*1e9,fit_tmp(fwn1:fwn2),Lstyle_tmp,'DisplayName','Fit')
            end

            % Plot the ratio of the derivatives
            plot(Dp_mean_bin(fwn1:fwn2).*1e9,abs(ddDp_ratio_Dp2molec_tmp(fwn1:fwn2)),...
                'LineStyle',Lstyle,'Color',colors_for_runs(nrun,:),'LineWidth',Lwidth+1,'DisplayName','\partial^2:\partial')

            for nth=1:length(fig_th)
                if nth < length(fig_th)
                    plot([Dp_limits_bin(1)*1e9 Dp_limits_bin(end)*1e9],[fig_th(nth) fig_th(nth)],':','Color','k','LineWidth',Lwidth,'HandleVisibility','off')
                else
                    plot([Dp_limits_bin(1)*1e9 Dp_limits_bin(end)*1e9],[fig_th(nth) fig_th(nth)],':','Color','k','LineWidth',Lwidth,'DisplayName',legend_for_th)
                end
            end

            xlim([Dp_limits_bin(fwn1)*0.5 Dp_limits_bin(fwn2+1)*2]*1e9)
            ylim_tmp=get(gca,'YLim');
            ylim([ylim_tmp(1)*0.5,ylim_tmp(2)*2])

            xlabel('{\itd}_{p,mob} (nm)')
            ylabel('{\itC} (cm^{-3}/m)')

            plot([Dp_limits_bin(nbin)*1e9 Dp_limits_bin(nbin)*1e9],ylim_tmp,'--','Color',colors_for_runs(nrun,:),'LineWidth',Lwidth,'HandleVisibility','off')
            plot([Dp_limits_bin(nbin+1)*1e9 Dp_limits_bin(nbin+1)*1e9],ylim_tmp,'--','Color',colors_for_runs(nrun,:),'LineWidth',Lwidth,'HandleVisibility','off')

            title_str_tmp=['Fit around bin ',legend_for_bins{nbin},' at ',sprintf('%.0f',time_vector(tind)),' min'];
            if tval == apptime_bin(nbin)
                title_str_tmp=[title_str_tmp,' ({\itt}_{app} of the bin)'];
            end
            title(title_str_tmp,'FontWeight','normal')
            legend('Location','best')
            
            if l_pause == 1
                pause
            end
            
        end
        
        % Save the value for each bin if no time-dependent metrics are assessed
        if ismember(time_tag,{'tapp','last'})
            ddDp_ratio_Dp2molec(nbin)=abs(ddDp_ratio_Dp2molec_tmp(nbin));
            r2(nbin)=gof.rsquare;
        end

        % Save the time-dependent derivative ratio
        ddDp_ratio_Dp2molec_vs_time(nt,nbin)=abs(ddDp_ratio_Dp2molec_tmp(nbin));
        r2_vs_time(nt,nbin)=gof.rsquare;

    end
    
end

if ibins_incl_min == ibins_incl_max
    str_tmp=[num2str(ibins_incl_min),' for all fits'];
else
    str_tmp=sprintf('%d-%d',ibins_incl_min,ibins_incl_max);
end
fprintf('\nNumber of data points per fit: %s\n',str_tmp)
if ibins_incl_min < ibins_incl_def
    fprintf('! Would be good with at least ~%d for all fits if the size resolution is high\n',ibins_incl_def)
end

% R2 below given threshold
l_r2_low=zeros(size(r2)); l_r2_low(r2 < r2_th)=1;
l_r2_vs_time_low=zeros(size(r2_vs_time)); l_r2_vs_time_low(r2_vs_time < r2_th)=1;

%%

% Analyse the trends in the derivatives

if ismember(time_tag,{'tapp','last'})
    
    % Derivatives at appearance times or for a single size distribution
    
    ddDp_ratio=ddDp_ratio_Dp2molec;
    
    figure(12)
    set(gcf,'Position',Fpos)
    set(gca,'LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize)
    set(gcf,'Color','white')
    set(gca,'YScale','log')
    box on
    hold on

    plot(Dp_mean_bin.*1e9,ddDp_ratio,'LineStyle',Lstyle,'Color',colors_for_runs(nrun,:),...
        'LineWidth',Lwidth+1,'DisplayName','\partial^2:\partial')

    % Mark the values for which the size distribution fit was worse
    ind_tmp=find(l_r2_low==1);
    plot(Dp_mean_bin(ind_tmp).*1e9,ddDp_ratio(ind_tmp),'LineStyle','None','Marker',Mstyle_point,'MarkerSize',Msize_point,...
        'Color',colors_for_runs(nrun,:),'LineWidth',Lwidth,'DisplayName',l_r2_low_str)
    
    for nth=1:length(fig_th)
        if nth < length(fig_th)
            plot([Dp_limits_bin(1)*1e9 Dp_limits_bin(end)*1e9],[fig_th(nth) fig_th(nth)],':','Color','k','LineWidth',Lwidth,'HandleVisibility','off')
        else
            plot([Dp_limits_bin(1)*1e9 Dp_limits_bin(end)*1e9],[fig_th(nth) fig_th(nth)],':','Color','k','LineWidth',Lwidth,'DisplayName',legend_for_th)
        end
    end
    
    xlim([Dp_limits_bin(nbin_start) Dp_limits_bin(nbin_end+1)]*1e9)

    xlabel('{\itd}_{p,mob} (nm)')
    ylabel('\partial^2:\partial')
    
    title_str_tmp='Value for each bin';
    if strcmp(time_tag,'tapp')
        title_str_tmp=[title_str_tmp,' at appearance time'];
    elseif strcmp(time_tag,'last') && length(time_vector) > 1
        title_str_tmp=[title_str_tmp,' at the end'];
    end
    title(title_str_tmp,'FontWeight','normal')
    legend('Location','best')

    % Find where the ratio of the 2nd and 1st derivatives of the distribution decreases below a given threshold
    % There's a placeholder for adding more variables, e.g. for finding where the relative difference between GR_app and GR_cond
    % decreases below a given threshold
    var_name={'d2_d'};
    nvar=length(var_name);                                  % Number of variables to be studied
    var_value=struct('data',{ddDp_ratio},'th',{small_th});  % Variables and thresholds

    % Use the following conditions to assess the decrease beyond the threshold
    cond_name={'FirstPass', 'LastPass', 'FirstPass_Ave', 'LastPass_Ave'};
    ncond_th=length(cond_name);                                         % Number of conditions to be checked
    cond_str={'var_data(nbin) < var_th',...                             % First occurence where the variable decreases below the threshold
        'max(var_data(nbin:nbin_max),[],''omitnan'') < var_th',...      % Point after which the variable doesn't increase above the threshold anymore
        'var_data_ave(nbin) < var_th',...                               % Same as above, but for averaged values
        'max(var_data_ave(nbin:nbin_max),[],''omitnan'') < var_th'};

    bin_th=nan(nvar,ncond_th);
    Dp_th=nan(nvar,ncond_th);

    for nv=1:nvar

        var_data=var_value(nv).data;
        var_th=var_value(nv).th;

        % See also averaged data (running mean over given number of points along diameter axis)
        var_data_ave=smooth(var_data,5);

        lfound=zeros(1,ncond_th);

        for nbin=1:ibins

            for ncond=find(lfound==0)
                eval(['cond_true = ',cond_str{ncond},';']);
                if cond_true
                    bin_th(nv,ncond)=nbin;
                    Dp_th(nv,ncond)=Dp_mean_bin(nbin);
                    lfound(ncond)=1;
                end
            end

            if all(lfound==1)
                break
            end

        end

        % Plot also the averaged quantities to see if the averaging makes sense
        figure(12)
        ind=find(var_data_ave>=0);
        plot(Dp_mean_bin(ind).*1e9,var_data_ave(ind),'LineStyle',Lstyle,'Color',[.64 .08 .18],'LineWidth',Lwidth,'DisplayName','Running mean')

    end

    % Print out a table showing the sizes corresponding to the thresholds
    table_th = array2table(Dp_th*1e9,'RowNames',var_name,'VariableNames',cond_name);
    fprintf('\n-----------------------------------------------------------------------------------------------\n\n')
    disp('Convergence size (nm)')
    fprintf(['Convergence criterion: ',num2str(small_th*100),'%%\n\n'])
    disp(table_th)
    
else
    
    % Assessment of the size and time ranges where GR is expected to be applicable
    
    % Shift time axis labels to be in the middle of the patches
    times_diff=diff(times_for_ddDp_ratio);
    times_shifted=times_for_ddDp_ratio(1:end-1)+times_diff/2;
    times_shifted=[times_shifted times_for_ddDp_ratio(end)+times_diff(end)/2];
    
    gr_cond_ok=zeros(size(ddDp_ratio_Dp2molec_vs_time));
    l_r2_points=[];
    
    for nt=1:nt_tot
        for nbin=1:ibins
            if ddDp_ratio_Dp2molec_vs_time(nt,nbin) < small_th
                gr_cond_ok(nt,nbin)=1;
            end
            if l_r2_vs_time_low(nt,nbin) == 1
                l_r2_points=[l_r2_points; times_shifted(nt), Dp_mean_bin(nbin)];
            end
        end
    end
    
    figure(1)
    clf(1)
    set(gcf,'Position',Fpos)
    set(gcf,'Color','white');
    box on
    hold on
    
    % Add extra elements to the data as pcolor / surf discards the last row and column
    xaxis_plot=[times_for_ddDp_ratio times_for_ddDp_ratio(end)+times_diff(end)];
    yaxis_plot=Dp_limits_bin*1e9;
    
    data_plot=gr_cond_ok';
    data_plot=[data_plot nan(size(data_plot,1),1); nan(1,size(data_plot,2)+1)];
    
    h_tmp=pcolor(xaxis_plot,yaxis_plot,data_plot);
    h_tmp.Annotation.LegendInformation.IconDisplayStyle='off';
    set(gca,'LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize)
    colormap([1 1 1; 0.8 0.9 0.5])
    caxis([0 1])
    shading flat
    
    title('Shaded area: Continuous regime, GR expected to be ok','FontWeight','normal')
    
    % Show d2:d as a function of size and time
    figure(2)
    clf(2)
    set(gcf,'Position',Fpos)
    set(gcf,'Color','white')
    box on
    hold on
    
    data_plot=log10(ddDp_ratio_Dp2molec_vs_time');
    data_plot=[data_plot nan(size(data_plot,1),1); nan(1,size(data_plot,2)+1)];
    
    h_tmp=pcolor(xaxis_plot,yaxis_plot,data_plot);
    h_tmp.Annotation.LegendInformation.IconDisplayStyle='off';
    set(gca,'LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize);
    colormap(jet)
    shading flat
    
    cb=colorbar('LineWidth',Lwidth,'FontWeight','normal','FontSize',Fsize);
    lim_tmp=[floor(log10(small_th)) 1];
    caxis(lim_tmp)
    tick_tmp=unique([floor(min(min(data_plot,[],'omitnan'))):ceil(max(max(data_plot,[],'omitnan'))) lim_tmp]);
    set(cb,'YTick',tick_tmp)
    set(cb,'YTickLabel',arrayfun(@(x) num2str(x,'10^{%d}'), tick_tmp,'UniformOutput',false))
    ylabel(cb,'\partial^2:\partial','FontWeight','normal','FontSize',Fsize);

    for nfig=1:2
        
        figure(nfig)
        
        if ~isempty(l_r2_points)
            
            plot(l_r2_points(:,1),l_r2_points(:,2)*1e9,'LineStyle','None','Marker',Mstyle_point,'MarkerSize',Msize_point,'Color','k',...
                'LineWidth',Lwidth,'DisplayName',l_r2_low_str)
            
            legend('Location','Best')
            
        end
        
        set(gca,'Layer','top')
        
        % Limit the number of tick labels for figure clarity
        n_ticks=10;
        if length(times_for_ddDp_ratio) > n_ticks
            int_tmp=ceil(length(times_for_ddDp_ratio)/n_ticks);
        else
            int_tmp=1;
        end
        
        set(gca,'XTick',times_for_ddDp_ratio(1:int_tmp:end)); xtickformat('%.0f')
        label_tmp=get(gca,'XTickLabel');
        set(gca,'XTick',times_shifted(1:int_tmp:end)); set(gca,'XTickLabel',label_tmp)
        
        xlabel('{\itt} (min)')
        ylabel('{\itd}_{p,mob} (nm)')
        xlim([xaxis_plot(1)-1 xaxis_plot(end)+1])
        ylim([Dp_limits_bin(nbin_start)*1e9 Dp_limits_bin(nbin_end+1)*1e9])
        
    end
    
end

drawnow

end
