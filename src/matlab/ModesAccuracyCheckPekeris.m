function  err = ModesAccuracyCheckPekeris(wnum,MP,freq)

omeg = 2*pi*freq;
H = MP(2,1);
cw = MP(2,2);
cb = MP(2,3);
rhow = MP(2,4);
rhob = MP(2,5);

kvw = sqrt( (omeg/cw)^2 - wnum.^2 );
kvb = sqrt(  wnum.^2 - (omeg/cb)^2 );


err = tan( kvw*H ) + rhob*kvw./( rhow*kvb );