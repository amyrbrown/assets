
%
% omm3.m
%
% Derives an expression for the steady-state of the Open
% Michaelis-Menten model using a substitution strategy that 
% results in a sublinear term, and illustrates the use of a 
% pseudospecies to resolve the non-linearity.
%
% Optional Arguments:
% doUnsub - if true, the py-substitution is reversed prior
%       to returning vbar. Default value is true.
% verbose - if true, each calculation in the substitution
%       and solution of the steady state equation are 
%       printed to the console. Default value is false.
% fid - a file id to which to print text when verbose=true. 
%       Returned by a previous call to fopen(). Default 
%       value is 1, the console. This argument is ignored 
%       if verbose=false.
% wid - width of the device, in characters, to which to
%       print. Default value is 80.
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

function [vbar,psi] = omm3( varargin )

  ip = inputParser;
  addParamValue(ip,'doUnsub',true,@islogical);
  addParamValue(ip,'verbose',false,@islogical);
  addParamValue(ip,'fid',1,@isnumeric);
  addParamValue(ip,'wid',80,@isnumeric);

  parse(ip,varargin{:});
  doUnsub = ip.Results.doUnsub;
  verbose = ip.Results.verbose;
  fid = ip.Results.fid;
  wid = ip.Results.wid;

  % Load the model description
  [S,v,k,x,xdot] = loadOMM();

  if( verbose )
    printvar(sprintf(strcat('\nThere are %d reactions',...
      ' and %d species in the OMM model.\n'),...
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

  % Add the pseudospecies x5_hat=1
  syms x5_hat positive real
  x = [ x; x5_hat ];
  v = subs(v,k(4),k(4)*x5_hat);
  
  if( verbose )
    printvar(sprintf(strcat('\nWe would like to',...
      ' define a map psi_p such that psi_p(k4) is in P.',...
      ' To do so we introduce a pseudospecies x5_hat=1',...
      ' and let v(4) = k4*x5_hat. This gives\n')),...
      'wid',wid,'fid',fid);
    printvar( v, 'vsym', 'v', 'fid', fid, 'wid', wid );
  end

  % Declare vectors of in/dependent variables
  p = sym( 'p', [5,1] ); y = sym( 'y', [5,1] );
  p = sym( p, 'positive' ); p = sym( p, 'real' );
  y = sym( y, 'positive' ); y = sym( y, 'real' );
  
  % Substitution strategy
  kx = [ k' x' ];
  py = [ y(1:3)',p(5),y(4),p(1:4)',y(5) ];
  
  % Execute pysub
  [vbar,ybar,psi_py,q] = pysub( S, v, kx, py, y,...
    'verbose', verbose, 'fid', fid, 'wid', wid );

  % Resolve the pseudospecies
  vbar = subs(vbar,psi_py(find(psi_py==x5_hat),2),1);
  ybar = subs(ybar,psi_py(find(psi_py==x5_hat),2),1);
  psi_py_ok = subs(psi_py,psi_py(find(psi_py==x5_hat),2),1);

  if( verbose )
    lhs = char(psi_py(find(psi_py==x5_hat),2));
    streq = strcat( lhs, '=1' );
    soln = strcat(lhs,'=',char(solve(streq,lhs)));
    printvar(sprintf(strcat('\nTo resolve the pseudospecies',...
      ' we require that psi_py(x5_hat)=1. In other words, %s.',...
      ' This gives\n'), soln),'wid',wid,'fid',fid);
    printvar( psi_py_ok, 'vtype', 'map', 'vsym', 'psi_py',...
      'fid', fid, 'wid', wid );
    printvar(sprintf('\nand\n'),'wid',wid,'fid',fid);
    printvar( vbar, 'vtype', 'vec', 'vsym', 'vbar',...
      'fid', fid, 'wid', wid );
  end
  
  % Reverse the substitution
  if( doUnsub ) 
    if( verbose )
      printvar(sprintf(strcat('\nWe may now proceed with',...
        ' the inverse substitution.')),'wid',wid,'fid',fid);
    end
    [vbar,psi_ss] = unsub( v,p,q,y,kx,py,ybar,psi_py_ok,...
      'verbose', verbose, 'fid', fid, 'wid', wid );
    psi = psi_ss;
  else
    psi = psi_py_ok;
  end
  
  % Verify steady state
  if( verbose )
    printvar(sprintf(strcat('\nVerify steady state.',...
      ' xdot = S * vbar, where\n')), 'fid', fid, 'wid', wid); 
    printvar( simplify(S*vbar), 'vtype', 'vec', 'vsym', 'xdot',...
      'fid', fid, 'wid', wid );
  end
