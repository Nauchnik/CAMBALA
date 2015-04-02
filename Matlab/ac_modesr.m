function [wnum, wmode, varargout] = ac_modesr(z,MediaParams,freq, varargin)

% Solution of acoustical spectral problem + Richardson extrapolation;

%% mode functions from the solution on a given grid

if nargout == 3
    if nargin == 4
        [wnum, wmode, dwmode] = ac_modes(z,MediaParams,freq, varargin{1});
    else
        [wnum, wmode, dwmode] = ac_modes(z,MediaParams,freq);
    end;
    varargout{1} = dwmode;
else
    if nargin == 4
        [wnum, wmode] = ac_modes(z,MediaParams,freq, varargin{1});
    else
        [wnum, wmode] = ac_modes(z,MediaParams,freq);
    end;
end;

%% increasing the accuracy for wavenumbers

% Solving spectral problem for several grids with different meshsize

Ngr = 3; % number of grids
nmod = length(wnum);

dz0 = z(2) - z(1);  % "old" meshsize
zc = z;             % z-mesh

% initializing arrays 

wnums(1:Ngr,1:nmod) = repmat(wnum,Ngr,1); % wnumbers for different meshsizes
dzs(1:Ngr) = dz0; % meshsizes


for ii = 1:Ngr
    
    % new meshsizes are chosen so that
    % they're close to the 1.5*dz, 2*dz, 2.5*dz, etc
    % but bottom depth is a multiple of meshsize
    
    
    dscDepths = MPFindDsc(MediaParams);
    
    
    if ~isempty(dscDepths)
        bDepth = dscDepths(1);
    else
        bDepth = z(end);
    end;
    
     %izH = fix( bDepth/(dz0*(0.5+ii/2) ) ); 
     izH = fix( bDepth/(dz0*(ii) ) ); 
     
        if izH > 0
            dzc = bDepth/izH;
            zc = 0:dzc:z(end);
        end;
    
    dzs(ii) = dzc;    
    [wnum, ~] = ac_modes(zc,MediaParams,freq, nmod);
    wnums(ii,1:nmod) = wnum(1:nmod);
    

end;

% Richardson extrapolation

rMatrix(1:Ngr,1:Ngr) = repmat(dzs(1:Ngr).',1,Ngr).^repmat( 2*(0:Ngr-1),Ngr,1);
rDat = rMatrix\(wnums.^2);
wnum = sqrt( rDat(1,1:nmod) );

AA = rMatrix^(-1);
disp('Richardson coeffs:');
disp( AA(1,:) );


