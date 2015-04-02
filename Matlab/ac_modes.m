function [wnum, wmode, varargout] = ac_modes(z,MediaParams,freq, varargin)

% In the eigenvectors matrix Vec (wmode) modes are written as colums, i.e.
% first mode is Vec(:,1) , the second one is Vec(:,2) and so on.
% freq  is frequency in Hertzs;
% y is depths vector, y(1)=0 - surface, y(ny) - bottom depth 
% where Neumann boundary condition u_z=0 is imposed

% syntax1: [wnum, wmode] = ac_modes(z,MediaParams,freq) - compute wave
% numbers, acoustic modes for given depth vector, MediaParams and frequency

% syntax2: [wnum, wmode, dwmode] = ac_modes(z,MediaParams,freq)
% compute derivatives of acoustic modes - dwmode - as well 

% syntax3: [wnum, wmode, dwmode] = ac_modes(z,MediaParams,freq,nmod)
% compute only first nmod wave numbers and modes, for nmod = 0 function
% computes only propagating modes

if isstruct(MediaParams)
    LayersData = MediaParams.LayersData;
else
    LayersData = MediaParams;
end;



omeg = 2*pi*freq;
nz = length(z);
dz = z(2) - z(1);
nDsc = length(LayersData(:,1));
ddepths(1:nDsc) = LayersData(1:nDsc,1);
cm(1:nDsc) = LayersData(1:nDsc,2);
cp(1:nDsc) = LayersData(1:nDsc,3);
dm(1:nDsc) = LayersData(1:nDsc,4);
dp(1:nDsc) = LayersData(1:nDsc,5);
nmod = nz;
%% constructing sound speed profile
c(1:nz)=cm(1);
ziParams(1:nDsc) = round(LayersData(1:nDsc,1)/dz)+1;


for ii=2:nDsc
    if ziParams(ii-1)~=ziParams(ii)
        c(ziParams(ii-1):ziParams(ii))=interp1([ddepths(ii-1), ddepths(ii)],...
            [cp(ii-1), cm(ii)],z(ziParams(ii-1):ziParams(ii)),'linear','extrap');
    else
        c(ziParams(ii)) = cp(ii);
    end;
end;

c(ziParams(nDsc):end)=LayersData(nDsc,3);

if isstruct(MediaParams)
    if isfield(MediaParams,'HydrologyData')
        if nDsc > 1
            if ziParams(2) > 1
                c(1:ziParams(2)) = interp1(MediaParams.HydrologyData(:,1),MediaParams.HydrologyData(:,2),z(1:ziParams(2)),'linear','extrap');
                cm(2) = c(ziParams(2));
                c(ziParams(2)) = cp(2);
            end;
        else
            c = interp1(MediaParams.HydrologyData(:,1),MediaParams.HydrologyData(:,2),z,'linear','extrap');
        end;
    end;
end;



%% constructing density profile
d(1:nz)=dm(1);
ziParams(1:nDsc) = round(LayersData(1:nDsc,1)/dz)+1;

for ii=2:nDsc
    if ziParams(ii-1)~=ziParams(ii)
        d(ziParams(ii-1):ziParams(ii))=interp1([ddepths(ii-1), ddepths(ii)],...
            [dp(ii-1), dm(ii)],z(ziParams(ii-1):ziParams(ii)),'linear','extrap');
    else
        d(ziParams(ii)) = dp(ii);
    end;
end;

d(ziParams(nDsc):end)=LayersData(nDsc,5);

%% assembling tridiagonal matrix

% treating continuity points

lower(1:nz) = 1/(dz^2);
upper(1:nz) = 1/(dz^2);
main(1:nz) = -2/(dz^2) + (omeg./c(1:nz)).^2;

% Neumann boundary condition at the bottom

lower(nz) = 4*d(nz)/((dz^2)*( d(nz-1) + d(nz) ));
main(nz) = - ( d(nz-1) + 3*d(nz) )/( (dz^2)*( d(nz) + d(nz-1) )) + (omeg./c(nz)).^2;

% Dirichlet pressure release at the top

lower(1) = 0;
lower(2) = 0;
upper(1) = 0;

%treating discontinuities of media params

stDsc = find(ziParams>=2,1,'first');

lower(ziParams(stDsc:nDsc)) = (2/(dz^2))*(dp(stDsc:nDsc)./(dp(stDsc:nDsc)+dm(stDsc:nDsc)));
main(ziParams(stDsc:nDsc)) = (omeg^2)*(dm(stDsc:nDsc)./((dp(stDsc:nDsc)+dm(stDsc:nDsc)).*cp(stDsc:nDsc).^2 ) +...
        dp(stDsc:nDsc)./((dp(stDsc:nDsc)+dm(stDsc:nDsc)).*cm(stDsc:nDsc).^2 ) ) - 2/(dz^2);
upper(ziParams(stDsc:nDsc)) = (2/(dz^2))*(dm(stDsc:nDsc)./(dp(stDsc:nDsc)+dm(stDsc:nDsc)));


dsm = spdiags([upper' main' lower'],-1:1,nz,nz )';



%% solution of the acoustic spectral problem

[Vec,Dia] = eig(full(dsm));

[wvect,indx] = sort(diag(Dia,0));
Vec = Vec(1:nz,indx);
wnumAll(1:nmod) = sqrt(flipud(wvect(nz-nmod+1:nz) ));

% nmod =
%       0 -- propagating modes;
%       -1 -- all modes;
%       n -- n modes
%       not provided -- trapped modes


if nargin>3 
    if varargin{1}==0
        nmod = find(imag(wnumAll)==0,1,'last');
    elseif varargin{1}>0
        nmod = varargin{1};
    end;
else
    nmod = find(real(wnumAll) > omeg/cp(end),1,'last');
end;


wnum(1:nmod) = sqrt(flipud(wvect(nz-nmod+1:nz) ));
wmode(1:nz,1:nmod) = Vec(1:nz,nz-nmod+1:nz);

wvectf(1:nmod) = wmode(10,1:nmod);
wvectfm = repmat(wvectf,nz,1);

wmode = fliplr(wmode.*sign(wvectfm));
mnorms(1:nmod) = dz*trapz((wmode.^2)./repmat(d',1,nmod));
wmode = wmode./repmat(mnorms.^(0.5),nz,1);

%% mode derivatives

if nargout>2
    if ziParams(2)-ziParams(1)~=1
        dwmode(1,1:nmod) = ( 4*wmode(2,1:nmod) - 3*wmode(1,1:nmod) - wmode(3,1:nmod) )/(2*dz);
    else
        dwmode(1,1:nmod) = (wmode(2,1:nmod) - wmode(1,1:nmod))/dz;
    end;
    
    
    
    dwmode(2:nz-1,1:nmod) = (wmode(3:nz,1:nmod) - wmode(1:nz-2,1:nmod))/(2*dz);
    dwmode(nz,1:nmod) = 0;

    iDsc = ziParams(2:end);
    dwmode(iDsc,1:nmod) = ( 4*wmode(iDsc+1,1:nmod) - 3*wmode(iDsc,1:nmod) - wmode(iDsc+2,1:nmod) )/(2*dz);

    
    varargout{1} = dwmode;
end;
