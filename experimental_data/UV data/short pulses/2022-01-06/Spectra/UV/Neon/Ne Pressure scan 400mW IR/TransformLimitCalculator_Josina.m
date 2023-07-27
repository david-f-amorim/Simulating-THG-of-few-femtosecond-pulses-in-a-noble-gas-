clear all;
close all;


%% REPLACE COMMA DELIMITER FOR POINT DELIMITER %%

% replaces all occurences of comma (",") with point (".") in a text-file.
% Note that the file is overwritten, which is the price for high speed.
% file_withcomma = uigetfile('*.*');
% comma2point_overwrite( file_withcomma )


%% IMPORT DATA %%

file=uigetfile('*.*');

raw = importdata(file);
spectrum = raw.data;

wavelength = spectrum(:,1);
intensity_w = spectrum(:,2)./max(spectrum(:,2));

%% PLOT RAW SPECTRUM %%

figure('name','Raw spectrum');
plot(wavelength,intensity_w);
ylim([0 1.1]);
xlim([min(wavelength) max(wavelength)]);
xlabel('Wavelength (nm)');
ylabel('Intensity (norm.)');
title('Raw spectrum');
pbaspect([1 1 1])

%% ASK TO SELECT SPECTRAL REGION %%

mydialog('Spectral region selection','Select the relevant spectral region','OK')

[x,~] = ginput(2);
w_min = min(x);
w_max = max(x);
clear x;

[~,w_min_i] = min(abs(wavelength-w_min));
[~,w_max_i] = min(abs(wavelength-w_max));

wavelength = wavelength(w_min_i:w_max_i);
intensity_w = intensity_w(w_min_i:w_max_i);
intensity_w = intensity_w-min(intensity_w);

close 'Raw spectrum'

figure('name','Transform limited calculator');
%set(gcf, 'Position', [100, 100, 600, 1000])
subplot(121);
area(wavelength,intensity_w./max(intensity_w),'Linewidth',2, 'EdgeColor',[0.273 0.316 0.972]', 'FaceColor',[0.719 0.734 0.977]);
% '#4DBEEE'
ylim([0 1.15]);
xlim([min(wavelength) max(wavelength)]);
xlabel('Wavelength (nm)');
ylabel('Intensity (norm.)');
title(['Spectrum']);
% title([file ' - Relevant spectrum']);
pbaspect([1 1 1])


%% CALCULATE CARRIER WAVELENGTH %%

carrier_wavelength = sum(wavelength.*intensity_w)/sum(intensity_w) %% in [nm]


h = 4.13566751e-15; % in [eV][s]
c = 299792458;      % in [m][s]-1
carrier_photon_energy = h*c/carrier_wavelength*1e9 %% in [eV]


%% PERFORM THE FOURIER TRANSFORM%%

wavelength_nm = wavelength*1e-9;
f = 3e8./wavelength_nm;                     % calculates the frequency axis

profile_w = sqrt(intensity_w);              % calculates the spectral profile

[t,profile_t]=DFT(0,f,profile_w,0,0);

temp=[-max(t):abs(t(2)-t(1)):max(t)]; 
t=temp; clear temp;                          % creates symmetric time axis

temp=zeros(1,2*length(profile_t)-1);
temp(1:length(profile_t))=fliplr(profile_t);
temp(length(profile_t):end)=profile_t;
profile_t=temp; clear temp;                 % creates symmetric profile_w axis
                                       

intensity_t = (abs(profile_t)).^2;
intensity_t = intensity_t./max(intensity_t);

set(gca,'FontSize',22)
set(gca,'fontname','times')

subplot(122);
plot(t.*1e15,intensity_t,'.','Markersize',10);
xlim([-20 20])
ylim([0 1.15])
xlabel('Time (fs)');
ylabel('Intensity (norm.)');
title('Time domain intensity');
pbaspect([1 1 1])
%% GAUSSIAN FIT %%

[parameters,fitness] = fit(t'.*1e15,intensity_t','gauss1');
t_fit = [-20:0.001:20];
a1 = parameters.a1;
b1 = parameters.b1;
c1 = parameters.c1;
fit = a1*exp(-((t_fit-b1)/c1).^2);
hold on
plot(t_fit,fit,'Linewidth',1.5, 'color', 'red')
sigma = parameters.c1/sqrt(2);
TL = 2.355*sigma;
set(gcf,'color','w');
legend('Fourier transform of the spectrum','Gaussian fit to the fourier transform')

% dim = [.6 .6 .3 .3];
dim = [0.75 0.7 0.05 0.1];
str = ['The TL duration is ',num2str(TL),' fs'];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'fontname','times')
set(gca,'FontSize',22)

mydialog('Transform Limit Duratoin',['The TL duration is ',num2str(TL),' fs'],'OK')

%% plot spectrum nice

% figure('name','Transform limited calculator');
% %set(gcf, 'Position', [100, 100, 600, 1000])
% plot(wavelength,intensity_w./max(intensity_w), linewidth=1.5);
% ylim([0 1.15]);
% xlim([min(wavelength) max(wavelength)]);
% xlabel('Wavelength (nm)');
% ylabel('Intensity (norm.)');
% % title([strrep(file(1:end-4), '_', ' '), ' - Relevant spectrum']);
% pbaspect([1 1 1])
% set(gcf,'color','w');

%% DEFINE INPUT DIALOG %%

function mydialog(name,text,buttontext)

 d = dialog('Position',[500 500 250 150],'Name',name);

    txt = uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 210 40],...
               'String',text);

    btn = uicontrol('Parent',d,...
               'Position',[85 20 70 25],...
               'String',buttontext,...
               'Callback','delete(gcf)');
  
end


function    comma2point_overwrite( filespec )
% replaces all occurences of comma (",") with point (".") in a text-file.
% Note that the file is overwritten, which is the price for high speed.
    file    = memmapfile( filespec, 'writable', true );
    comma   = uint8(',');
    point   = uint8('.');
    file.Data( transpose( file.Data==comma) ) = point;
    %%%delete(file)
end
