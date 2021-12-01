function [ fea, out ] = ex_linearelasticity3( varargin )
%EX_LINEARELASTICITY3 Parametric study for the bracket deflection model.
%
%   [ FEA, OUT ] = EX_LINEARELASTICITY3( VARARGIN ) Example to conduct parametric
%   study of maximum stresses for the bracket deflection model.
%
%       Input       Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       thickness   [0.035 0.03 0.025 0.02]   Plate thickness
%       loads       [1e4 2e4 1e5 1e6]         Applied loads.
%                                                                                         .
%       Output      Value/(Size)           Description
%       -----------------------------------------------------------------------------------
%       fea         struct                 Problem definition struct
%       out         struct                 Output struct

% Copyright 2013-2021 Precise Simulation, Ltd.

thickness = [ 0.035 0.03 0.025 0.02 ];
loads     = [ 1e4 2e4 1e5 1e6 ];

%% Starting new model.
fea.sdim = { 'x', 'y', 'z' };
fea = addphys( fea, @linearelasticity, { 'u', 'v', 'w' } );

for i=1:length(thickness)

  %% Geometry operations.
  fea.geom = struct;
  gobj = gobj_block( 0, 0.02, 0, 0.2, 0, 0.2, 'B1' );
  fea = geom_add_gobj( fea, gobj );
  gobj = gobj_block( 0, 0.2, 0, 0.2, 0.1-thickness(i)/2, 0.1+thickness(i)/2, 'B1' );
  fea = geom_add_gobj( fea, gobj );
  gobj = gobj_cylinder( [ 0.1, 0.1, 0.08 ], 0.08, 0.04, 3, 'C1' );
  fea = geom_add_gobj( fea, gobj );
  fea = geom_apply_formula( fea, 'B1+B2-C1' );

  %% Grid operations.
  fea.grid = gridgen( fea, 'hmax', thickness(i)/2 );

  for j=1:length(loads)

    %% Equation settings.
    fea.phys.el.dvar = { 'u', 'v', 'w' };
    fea.phys.el.sfun = { 'sflag1', 'sflag1', 'sflag1' };
    fea.phys.el.eqn.coef = { 'nu_el', 'nu', 'Poissons ratio', { '0.3' };
                             'E_el', 'E', 'Modulus of elasticity', { '200e9' };
                             'Fx_el', 'F_x', 'Body load in x-direction', { '0' };
                             'Fy_el', 'F_y', 'Body load in y-direction', { '0' };
                             'Fz_el', 'F_z', 'Body load in z-direction', { '0' };
                             'u0_el', 'u0', 'Initial condition for u', { '0' };
                             'v0_el', 'v0', 'Initial condition for v', { '0' };
                             'w0_el', 'w0', 'Initial condition for w', { '0' } };
    fea.phys.el.eqn.seqn = { 'u'' - ((E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*ux_x + (E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy_x + (E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz_x + (E_el/(1+nu_el)/2)*(uy_y+vx_y) + (E_el/(1+nu_el)/2)*(uz_z+wx_z)) = Fx_el', 'v'' - ((E_el/(1+nu_el)/2)*(uy_x+vx_x) + (E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux_y + (E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*vy_y + (E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz_y + (E_el/(1+nu_el)/2)*(vz_z+wy_z)) = Fy_el', 'w'' - ((E_el/(1+nu_el)/2)*(uz_x+wx_x) + (E_el/(1+nu_el)/2)*(vz_y+wy_y) + (E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux_z + (E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy_z + (E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*wz_z) = Fz_el' };
    fea.phys.el.eqn.vars = { 'von Mieses stress', 'sqrt(((E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz)^2+((E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz)^2+((E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*wz)^2-((E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz)*((E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz)-((E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz)*((E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*wz)-((E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz)*((E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*wz)+3*((E_el/(1+nu_el)/2)*(uy+vx))^2+3*((E_el/(1+nu_el)/2)*(vz+wy))^2+3*((E_el/(1+nu_el)/2)*(uz+wx))^2)';
                             'x-displacement', 'u';
                             'y-displacement', 'v';
                             'z-displacement', 'w';
                             'Stress, x-component', '(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz';
                             'Stress, y-component', '(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*vy+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*wz';
                             'Stress, z-component', '(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*ux+(E_el/((1+nu_el)*(1-2*nu_el))*nu_el)*vy+(E_el/((1+nu_el)*(1-2*nu_el))*(1-nu_el))*wz';
                             'Stress, xy-component', '(E_el/(1+nu_el)/2)*(uy+vx)';
                             'Stress, yz-component', '(E_el/(1+nu_el)/2)*(vz+wy)';
                             'Stress, zx-component', '(E_el/(1+nu_el)/2)*(uz+wx)';
                             'Displacement field', { 'u', 'v', 'w' } };

    %% Constants and expressions.
    fea.expr = { 'load', loads(j) };

    %% Boundary settings.
    n_bdr = max(fea.grid.b(3,:));
    ind_fixed = findbdr( fea, 'x<0.005' );
    ind_load  = findbdr( fea, 'x>0.19'  );
    fea.phys.el.bdr.sel = ones(1,n_bdr);

    % Fix right boundary (set zero Dirichlet BCs).
    n_bdr  = max(fea.grid.b(3,:));         % Number of boundaries.
    bctype = num2cell( zeros(3,n_bdr) );   % First set homogenous Neumann BCs everywhere.
    [bctype{:,ind_fixed}] = deal( 1 );     % Set Dirchlet BCs for right boundary.
    fea.phys.el.bdr.coef{1,5} = bctype;

    % Apply negative z-load to left outer boundary.
    bccoef = num2cell( zeros(3,n_bdr) );
    [bccoef{3,ind_load}] = deal( '-load' );
    fea.phys.el.bdr.coef{1,end} = bccoef;

    %% Solver call.
    fea = parsephys( fea );
    fea = parseprob( fea );

    fea.sol.u = solvestat( fea );

    %% Postprocessing.
    s_vm = fea.phys.el.eqn.vars{1,2};
    vm_stress = evalexpr( s_vm, fea.grid.p, fea );
    max_vm_stress(i,j) = max(vm_stress);

  end
end


%% Visualization.
[t,l] = meshgrid(thickness,loads);
surf( t, l, max_vm_stress, log(max_vm_stress) )
xlabel('Thickness')
ylabel('Load')
zlabel('Maximum stress')
view(45,30)


out = max_vm_stress;
if ( nargout==0 )
  clear fea out
end
