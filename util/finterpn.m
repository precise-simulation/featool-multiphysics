%FINTERPN Linear N-D interpolation with data from file.
%
%   [ Vq ] = FINTERPN( FILENAME(X1, X2, ..., Xn, V), Xq1, Xq2, ..., Xqn)
%   Interpolate data from FILENAME in query points Xq1, Xq2, ..., Xqn.
%
%   Interpolation sample point coordinates (X1, X2, ..., Xn) and
%   corresponding function values (V) must be provided by FILENAME in
%   either of the following formats
%
%       .csv - ASCII text file with columns of comma separated values
%       .txt - ASCII text file with columns of space separated values
%       .xls/xslx - Microsoft Excel table format
%
%   The data in FILENAME must be stored as column vectors X1, X2, ..., Xn, and
%   and V consistent to be used with NDGRID, for example of 3D interpolation
%
%       x1, y1, z1, v1
%       x2, y1, z1, v2
%       x1, y2, z1, v3
%       x2, y2, z1, v4
%       x1, y3, z1, v5
%       x2, y3, z1, v6
%             ...
%       x1, y1, zN, vM-5
%             ...
%       x2, y3, zN, vM
%
%   For the .mat file format, must contain a variable named "data"
%
%       .mat - MATLAB file with a variable named "data"
%
%   where "data" can be an array as defined above, or a cell array
%   containing the interpolation points and function values.
%
%       data = {x, y, z, v};  % x/y/z/v size[nx,ny,nz]
%
%   Examples:
%
%       1) Command line use (1D interpolation)
%
%       x = linspace(0,1,10)';  % Sample point coordinates
%       v = (1:10)';  % Values at sample points
%
%       basename = 'test_data';
%       data = {x, v};  % Cell data format used by mat file
%       save([basename,'.mat'], '-mat', 'data');
%
%       data = [x, v];  % Array data format used by csv/txt file
%       csvwrite([basename,'.csv'], data);
%       save([basename,'.txt'], '-ascii', 'data');
%
%       results = finterpn([basename,'.csv'], [0.2, 0.5])
%       results = finterpn([basename,'.mat'], [0.2, 0.5])
%       results = finterpn([basename,'.txt'], [0.2, 0.5])
%
%       2) Equation coefficient (2D interpolation)
%
%       x = linspace(0,1,10);
%       y = linspace(0,1,20);
%       [xx,yy] = ndgrid(x,y);
%       v = rand(size(xx));
%
%       basename = 'test_data';
%       data = [xx(:), yy(:), v(:)];
%       csvwrite([basename,'.csv'], data);
%
%       fea.sdim = {'x', 'y'};
%       fea.geom.objects = {gobj_rectangle};
%       fea.grid = rectgrid();
%       fea = addphys(fea,@poisson);
%       source_term = ['1+2*finterpn(''',basename,'.csv'',x,y)'];
%       fea.phys.poi.eqn.coef{3,4} = { source_term };
%
%       fea = parsephys(fea);
%       fea = parseprob(fea);
%
%       fea.sol.u = solvestat(fea);
%
%       subplot(1,2,1)
%       postplot(fea, 'surfexpr', 'u')
%       subplot(1,2,2)
%       postplot(fea, 'surfexpr', source_term)
%
%   See also CSVREAD, XLSREAD, LOAD, NDGRID, INTERPN

% Copyright 2013-2025 Precise Simulation, Ltd.
