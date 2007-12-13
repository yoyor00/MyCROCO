function x=readattribute(url)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the attribute of an OPENDAP dataset
%
%  Retry (100 times) in case of network failure.
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%argdods='-A -e -v';
%argdods='-A -e +v';
argdods='-A +v';
%
nmax=100;
%
%
x=[];
ntry=0;
while isempty(x) 
  if ntry>nmax
    error(['READATTRIBUTE: repeated failures after ',num2str(nmax),' queries'])
  end
  ntry=ntry+1;
  try
    x=loaddap(argdods,url);
  catch
    x=[];
    disp(['READATTRIBUTE: did not work at ',num2str(ntry),' try: lets try again.'])
  end
end
%
return
