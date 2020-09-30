function Hydrus(root, root_bil, nC,ir,y,sy)

if ~isdeployed
    clear all;
    root = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
    root_bil = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global'; %For grid calculations
    root_out = 'W:\share\Hydrus_DOWNQ'; %For grid output
    nC='AFR';   % Which continent
    ir='ISIMIP2B';   % Which RCP level
    y=3;    % Which R-file
    sy=2050;% Which analysis year (2050 or 2100)
end

% CONT={'AFR','ASIA','AUS','CAM','EUR','NAM','SAM'};

% RCP_name{1} = 'hist';
% RCP_name{2} = 'RCP26';
% RCP_name{3} = 'RCP45';
% RCP_name{4} = 'RCP60';
% RCP_name{5} = 'RCP85';
% RCP_name{6} = 'ISMIP2B';

if isdeployed
    %     nC = str2num(nC);
    %     ir = str2num(ir);
    y =  str2num(y);
    sy = str2num(sy);
end

if strcmp('hist',ir)
    root_data = sprintf('%s\\data\\%s\\ISIMIP\\R\\%s',root_bil,nC,ir);
elseif strcmp('ISIMIP2B',ir)
    root_data = sprintf('%s\\data\\%s\\ISIMIP2B\\R',root_bil,nC);
else
    root_data = sprintf('%s\\data\\%s\\ISIMIP\\R\\%s_%d',root_bil,nC,ir,sy);
end

%% Extract file list
fnames = dir(root_data);

for j=1:numel(fnames)
    R_names{j} = sprintf('%s\\%s',root_data,(fnames(j).name));
end

%% load ACC and DIR file
disp('Loading ACC');
fname = sprintf('%s\\data\\%s\\acc.mat', root_bil, nC);
load(fname);

disp('Loading fdir');
fname = sprintf('%s\\data\\%s\\fdir.mat', root_bil, nC);
load(fname);

disp('Loading adir');
fname = sprintf('%s\\data\\%s\\adir.mat', root_bil, nC);
load(fname);

%% Start script


% for y=3:numel(R_names)

clearvars -except root root_bil root_out nC ir sy RCP_name res CONT nc root_data fnames R_names acc fdir adir y

fprintf('%s\n',fnames(y).name(1:(end-4)))

%% Load runoff
disp('Load runoff');
load(R_names{y})

%% Create sorting index
disp('Sorting');
[sortacc, order] = sort(acc(:), 'ascend');
%[~, order] = sort(Z(:), 'descend');

%% compute discharge by downstream routing
disp('Start routing')
defdirs

[nr, nc] = size(fdir);
sz = size(fdir);

Q = zeros(nr,nc, 'single');
N = nr*nc;
pold=0;
for k = 1:N
    i = order(k);
    r = rem(i-1, nr)+1; % Much faster than [r,c] = ind2sub(sz, i);
    c = fix((i-1)/nr)+1;
    q = R(r,c); % local discharge
    for d = 1:8;
        rj = r+drow(d);
        cj = c+dcol(d);
        if rj<1,  continue, end;
        if rj>nr, continue, end;
        if cj<1,  continue, end;
        if cj>nc, continue, end;
        if adir(rj,cj) ~= d, continue, end;
        %assert(Q(rj,cj)>0);
        q = q + Q(rj,cj);
    end;
    Q(r,c) = q;
    p = fix(1000*k/N);
    if p>pold
        fprintf('%d%% Continent %s File %s.\n',p,nC,fnames(y).name(1:(end-4)));
        pold=p;
    end;
end

%% For local runs

if ~isdeployed
    if strcmp('ISIMIP2B',ir)
        matpath = fullfile(root, sprintf('data\\%s\\ISIMIP2B\\Q', nC));
        if ~isdir(matpath)
            mkdir(matpath);
        end
        disp('Saving')
        matfile = fullfile(root, sprintf('data\\%s\\ISIMIP2B\\Q\\%s_%s_Q.mat',nC,fnames(y).name(1:(end-4))));
        save(matfile,'Q');
    else
        matpath = fullfile(root, sprintf('data\\%s\\ISIMIP\\Q\\%s', nC,ir));
        if ~isdir(matpath)
            mkdir(matpath);
        end
        disp('Saving')
        matfile = fullfile(root, sprintf('data\\%s\\ISIMIP\\Q\\%s\\%s_%s_Q.mat',nC,ir,fnames(y).name(1:(end-4))));
        save(matfile,'Q');
    end
end

%% For grid runs
if isdeployed
    root_out = 'W:\share\Hydrus_DOWNQ'; %For grid output
    matpath = fullfile(root_out, sprintf('%s', nC));
    if ~isdir(matpath)
        mkdir(matpath);
    end
    
    disp('Saving')
    if strcmp('hist',ir)
        matfile = fullfile(root_out, sprintf('%s\\Q_%s_%s.mat', nC,ir,fnames(y).name(1:(end-4))));
    elseif strcmp('ISIMIP2B',ir)
        matfile = fullfile(root_out, sprintf('%s\\Q_%s.mat', nC,fnames(y).name(1:(end-4))));
    else
        if sy==2100
            matfile = fullfile(root_out, sprintf('%s\\Q_%s_%d_%s.mat', nC,ir,sy,fnames(y).name(1:(end-4))));
        else
            matfile = fullfile(root_out, sprintf('%s\\Q_%s_%d_%s.mat', nC,ir,sy,fnames(y).name(1:(end-9))));
        end
    end
    save(matfile,'-v7.3','Q');
end

% end
%
%% show it.
% acc=single(acc);
% figure(1);clf;imagesc(log(acc));axis image;colormap(flipud(gray));
% [rmis,cmis]=find(acc==max(acc(:)));
% figure(1);hold on
% plot(cmis,rmis,'.r','markersize',20)
% hold off
% %%
% figure(2);clf;imagesc(log(Q));axis image; colormap(flipud(gray));
% [rmax,cmax] = find(Q==max(Q(:)));
% figure(2);hold on
% plot(cmax,rmax,'.b','markersize',20)
% plot(cmis,rmis,'.r','markersize',20)
% hold off
%
% %%
% figure(3);clf;imagesc(fdir);axis image; colormap(gray);
