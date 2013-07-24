
%
% fum1.m
%
% Solves for the steady-state enzyme concentrations in the
% fumarase model using both King-Altman and py-substitution,
% then demonstrates their equivalence. Note that fum1() does
% NOT accept a value for the doUnsub flag, since reversing
% the substitution is required to show equivalence with the
% King-Altman solution.
%
% Optional Arguments:
% verbose - if true, each calculation in the substitution and
%       solution of the steady state equation are printed to the 
%       console. Default value is false.
% fid - a file id to which to print text when verbose=true. Returned
%       by a previous call to fopen(). Default value is 1, the console. 
%       This argument is ignored if verbose=false.
% wid - width of the device, in characters, to which to print.
%       Default value is 80.
%
% Value:
% vbar - the steady state vector of reaction velocities 
%        satisfying S * vbar = 0. Elements in vbar are 
%        either in terms of p and q (if doUnsub is false), 
%        or k and x (if doUnsub is true).
% psi  - the composite mapping function from the set 
%        union(K,X) to union(P,Q,R[P,Q]) if doUnsub is 
%        false, or to union(K_pq,X_pq) if doUnsub is true, 
%        satisfying psi(v) = vbar.

function [vbar,psi] = fum1( varargin )

  ip = inputParser;
  addParamValue(ip,'verbose',false,@islogical);
  addParamValue(ip,'fid',1,@isnumeric);
  addParamValue(ip,'wid',80,@isnumeric);

  parse(ip,varargin{:});
  verbose = ip.Results.verbose;
  fid = ip.Results.fid;
  wid = ip.Results.wid;

  % Load the model description
  [S,v,k,x,xdot] = loadFum();

  if( verbose )
    printvar(sprintf(strcat('\nThere are %d reactions',...
      ' and %d species in the fumarase model.\n'),...
      length(k), length(x)), 'fid', fid, 'wid', wid );
    printvar(sprintf(strcat('The stoichiometric matrix is S =\n')),...
      'fid', fid, 'wid', wid);
    printvar( S, 'vsym', 'S', 'fid', fid, 'wid', wid );
    printvar(sprintf(strcat('\nThe vector of reaction',...
      ' velocities is v, where\n')),'wid',wid,'fid',fid);
    printvar( v, 'vsym', 'v', 'fid', fid, 'wid', wid );
    printvar(sprintf(strcat('\nThe vector of mass balance',...
      ' equations is xdot = S*v, where\n')),'wid',wid,'fid',fid);
    printvar( xdot, 'vsym', 'xdot', 'fid', fid, 'wid', wid );
  end
  
  % Assume substrates are constant
  S5 = S(1:5,:);
  xp = x(1:5);
  if( verbose )
    printvar(sprintf(strcat('\nTo compare py-substitution',...
      ' with King-Altman, we assume that the substrate',...
      ' abundances x6...x9 are constant.',...
      ' This gives S_5 = S(1:5,:) =\n')),'fid',fid,'wid',wid);
    printvar( S5, 'vsym', 'S', 'fid', fid, 'wid', wid );
    printvar(sprintf('\nand x'' =\n'),'fid',fid,'wid',wid);
    printvar( xp, 'vsym', 'x''', 'fid', fid, 'wid', wid );
  end

  
  %
  % First we'll solve by King-Altman
  %

  % Alternative 2D indexing for linear models
  syms kp13 kp34 kp15 kp54 kp42 kp21 positive real;
  syms kp31 kp43 kp51 kp45 kp24 kp12 positive real;
  kp = [ kp13,kp34,kp15,kp54,kp42,kp21,...
    kp31,kp43,kp51,kp45,kp24,kp12 ]';

  if( verbose )
    fprintf( fid, '\n\n%%' );
    fprintf( fid, '\n%% Solve S_5 * v = 0 by King-Altman' );
    fprintf( fid, '\n%%\n\n' );
    fprintf(fid,sprintf(strcat( 'First substitute in the',...
      ' following transition rate constants.\n\n')));
    ksub = subs(v,x(1:5),ones(5,1));
    printvar( [kp ksub], 'vtype', 'map', 'strop', '=',...
      'fid', fid, 'wid', wid );
  end

  % The Matlab subs command is not very smart
  % re: symbolic substitution. Substitute out rate
  % constants for transition rate constants then
  % set all "substrate" species x6...x9 equal to 1.
  vp = subs( v, k, kp );
  vp = subs( vp, x(6:9), ones(4,1) );

  % The resulting transition matrix K
  K = S5 * jacobian( vp, xp );
  if( verbose )
    printvar(sprintf(strcat('\nNow S_5 * v = K * x'',',...
      ' where K = S_5 * jacobian(v,x'') =\n')),...
      'fid',fid,'wid',wid);
    printvar( K, 'vsym', 'K', 'fid', fid, 'wid', wid );
  end

  % Calculate the minors
  M1 = det( K([2,3,4,5],[2,3,4,5]) );
  M2 = det( K([1,3,4,5],[1,3,4,5]) );
  M3 = det( K([1,2,4,5],[1,2,4,5]) );
  M4 = det( K([1,2,3,5],[1,2,3,5]) );
  M5 = det( K([1,2,3,4],[1,2,3,4]) );
  Mtot = simplify(M1+M2+M3+M4+M5);

  % Calculate steady state enzyme ratios
  xbar_ka = simplify( [ M1, M2, M3, M4, M5 ]' / Mtot );

  % Substitute out transition rate constants
  xbar_ka = subs( xbar_ka, kp, [ ...
    k(1)*x(6), k(2)*x(7), k(3)*x(7), k(4)*x(6), k(5)*x(8), k(6),...
    k(7), k(8), k(9), k(10), k(11), k(12)*x(9) ] );

  if( verbose )
    printvar(sprintf(strcat(...
      'Solve K * x'' = 0 by the King-Altman method.',...
      ' The resulting steady state expression for each',...
      ' enzyme i has the form N_ka(i)/D_ka, where\n'...
    )), 'fid',fid,'wid',wid);
    [N_ka,D_ka] = numden(xbar_ka);
    printvar( N_ka, 'vsym', 'N_ka', 'fid', fid, 'wid', wid, ...
      'dbspc', true );
    printvar(sprintf('and\n'),'fid',fid,'wid',wid);
    printvar( D_ka(1), 'vsym', 'D_ka', 'fid', fid, 'wid', wid );
  end


  %
  % Now we'll solve by py-substitution
  %

  % Declare vectors of in/dependent variables
  p = sym( 'p', [16,1] ); y = sym( 'y', [5,1] );
  p = sym( p, 'positive' ); p = sym( p, 'real' );
  y = sym( y, 'positive' ); y = sym( y, 'real' );
  
  % Substitution strategy
  kx = [ k' x' ];
  py = [ p(1:12)' y' p(13:16)' ];
  if( verbose )
    fprintf( fid, '\n\n%%' );
    fprintf( fid, '\n%% Solve S_5 * v = 0 by py-substitution' );
    fprintf( fid, '\n%%\n' );
  end

  % Execute pysub
  [vbar,ybar,psi_py,q] = pysub( S5, v, kx, py, y,...
    'verbose', verbose, 'fid', fid, 'wid', wid );

  % Reverse the substitution
  [vbar,psi_ss] = unsub( v,p,q,y,kx,py,ybar,psi_py, ...
    'verbose', verbose, 'fid', fid, 'wid', wid );
  psi = psi_ss;


  % Impose sum(x1...x5) = 1.
  xbar_py = psi(13:17,2);
  esum = simplify(sum(xbar_py));
  x5 = solve(sprintf('%s = 1',char(esum)),'x5');
  psi = subs( psi, x(5), x5 );
  vbar = subs( vbar, x(5), x5 );
  xbar_py = simplify(subs(xbar_py, x(5), x5 ));

  if( verbose )
    printvar(sprintf(strcat('\n\nThe stoichiometric matrix S5',...
      ' has a rank-deficiency of 1. In King-Altman, this',...
      ' free variable is x_tot. Unless specified, x_tot',...
      ' has an implicit value of 1:\n\n',...
      '>> simplify( sum(xbar_ka) ) =\n\n',...
      '  %s\n'), char(simplify(sum(xbar_ka)))),'fid',fid,'wid',wid );

    printvar(sprintf(strcat('\nIn py-substitution, the free',...
      ' variable is x5. To equate the two solutions, we',...
      ' require that sum(x1...x5) = 1.\n\n',...
      '>> x5 = solve( sprintf(''%%s = 1'',',...
      ' char(sum(xbar_py))), ''x5'' );\n')),'fid',fid,'wid',wid);

    printvar(sprintf(strcat(...
      '\nThe resulting steady state expression for each',...
      ' enzyme i has the form N_py(i)/D_py, where\n'...
    )), 'fid',fid,'wid',wid);
    [N_py,D_py] = numden(xbar_py);
    printvar( N_py, 'vsym', 'N_py', 'fid', fid, 'wid', wid, ...
      'dbspc', true );
    printvar(sprintf('and\n'),'fid',fid,'wid',wid);
    printvar( D_py(1), 'vsym', 'D_py', 'fid', fid, 'wid', wid );

    printvar(sprintf('\nVerify the two solutions are equal.\n'),...
      'fid',fid,'wid',wid);
    printvar(sprintf('>> delta_x = simplify(xbar_ka-xbar_py)\n'),...
      'fid',fid,'wid',wid);
    printvar( simplify(xbar_ka-xbar_py), 'vtype', 'vec', 'vsym', 'delta_x',...
      'fid', fid, 'wid', wid );
  end
  
