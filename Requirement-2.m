function [V, T] = Koch_snowflake_3D(nb_it, printable_ready, option_display)
%% Koch_snowflake_3D : function to compute, display, and save 
% a 3D Koch snowflake at any iteration / depth level.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
%
% Syntax
%
% Koch_snowflake_3D(nb_it);
% Koch_snowflake_3D(nb_it, printable_ready);
% Koch_snowflake_3D(nb_it, printable_ready, option_display);
% [V, T] = Koch_snowflake_3D(nb_it, printable_ready, option_display);
%
%
% Description
%
% Koch_snowflake_3D(nb_it) computes and displays nb_it
% 3D_Koch snowflake included in the unitary sphere.
%
% Koch_snowflake_3D(nb_it, printable_ready) prevents from
% creating non manifold edges when printable_ready is set to *true /
% logical 1, and remove duplicated vertices and faces when it is set to
% false / logical 0. In this latter case, the model is lighter (less
% vertices, less faces), but at the cost of non manifoldness.
%
% Koch_snowflake_3D(nb_it, printable_ready, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V, T] = Koch_snowflake_3D(nb_it, printable_ready, option_display) saves the resulting
% vertex coordinates in the array V, and the triangulation in the array T.
%
%
% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - printable_ready : either logical, true/*false or numeric 1/*0.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1
%
% Compute and display the 3D Koch snowflake at iteration 1, with minimum vertex and face numbers
% 
% Koch_snowflake_3D(1);
%
%
% Example #2
%
% Compute and save the 3D Koch snowflake at iteration 2, 3D printable ready, no display
%
% [V,T] = Koch_snowflake_3D(2,true,false);


%% Inputs parsing
assert(nargin < 4,'Too many input arguments.');

if ~nargin
    nb_it = 3;
    printable_ready = false;
    option_display = true;
elseif nargin > 0
    assert(isnumeric(nb_it) && nb_it == floor(nb_it) && nb_it >= 0,'nb_it parameter value must be numeric positive or null integer.');
    if nargin > 1
        assert(islogical(printable_ready) || isnumeric(printable_ready),'printable_ready parameter type must be either logical or numeric.');
        if nargin > 2
            assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
        else
            option_display = true;
        end
    else
        printable_ready = false;
        option_display = true;
    end
end

warning('on');
if option_display && nb_it > 3
    warning('%s triangles to display ! Make sure your graphic card has enough memory.',num2str(24*12^nb_it))    
end
warning('off');


%% Body
% Summits of original tetrahedron (living in S(O,1))
V1 = [2*sqrt(2)/3 0 -1/3];
V2 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];
V4 = [0 0 1];

Tetra = tetra(V1,V2,V3,V4);
Tetra_array = Tetra;

Triangle = triangle(V1,V2,V3);
Triangle_array = Triangle;

new_Triangle_array = repmat(Triangle_array,[1 1 6]);

% Build stellated octahedron and split into triangles -24 in total- and return the triangles structure array
T_current = Tetra_array;

[V_new,F_new] = split_tetra(T_current);

% Create new tetrahedrons
for m = 1:size(F_new,1)
    
    new_triangle = triangle(V_new(F_new(m,1),:),V_new(F_new(m,2),:),V_new(F_new(m,3),:));    
    new_Triangle_array(:,:,m) = new_triangle;
    
end

Triangle_array = new_Triangle_array; % 33 = 3 x 11 triangles

% Loop on nb_it
p = 0;

while p ~= nb_it
    
    new_Triangle_array = repmat(Triangle_array,[1 1 11]);       
    
    for j = 1:size(Triangle_array,3)
        
        T_current = Triangle_array(:,:,j);        
        [V_new,F_new] = split_triangle(T_current);
        
        for m = 1:size(F_new,1)
            
            new_triangle = triangle(V_new(F_new(m,1),:),V_new(F_new(m,2),:),V_new(F_new(m,3),:));            
            new_Triangle_array(:,:,11*(j-1)+m) = new_triangle;
            
        end                
        
    end
    
    Triangle_array = new_Triangle_array;    
    
    p = p+1;
    
end

% Triangles set concatenation
[V,T] = catriangles(Triangle_array);

if ~printable_ready
    
    % Remove duplicated vertices
    [V,T] = remove_duplicated_vertices(V,T);
    
    % Remove duplicated triangles
    T = unique(sort(T,2),'rows','stable');
    
end

% Display
if option_display
    
    cmap = [0 1 1];
    disp_koch_snowflake(V,T,cmap);
    
end


end % Koch_snowflake_3D


%% triangle structure computation subfunction
function [T] = triangle(V1, V2, V3)
%
% V1, V2, V3 : line vectors of triangle three vertices coordinates
%
% Vn = [Vxn Vyn Vzn]
%
% G = [Gx Gy Gz]

F = [1 2 3]; % summit / top first, then trigonometric order 

% triangle barycentre
G = mean([V1; V2; V3], 1);

T = struct('vertex', [V1; V2; V3], ...
           'centre', G, ...
           'facet', F);
       
end


%% tetrahedron structure computation subfunction
function [T] = tetra(V1, V2, V3, V4)
%
% V1, V2, V3, V4 : line vectors of tetrahedron four vertices coordinates
%
% Vn = [Vxn Vyn Vzn]
%
% G = [Gx Gy Gz]

F1 = [1 2 3];
F2 = [1 3 4];
F3 = [1 4 2];
F4 = [2 4 3];

% tetrahedron barycentre
G = mean([V1; V2; V3; V4], 1);

T = struct('vertex', [V1; V2; V3; V4], ...
           'centre', G, ...
           'facet', [F1; F2; F3; F4]);
       
end


%% Split triangle subfunction
function [V_new, F_new] = split_triangle(T)
%
% Input
%
% T : triangle
%
% Outputs
%
% V_new : The 11 vertices -coordinates- (8 new + 3 old)
% F_new : 11 indices line triplets of the 11 newly created facets


one_third = 1/3;
two_third = 2/3;

% Compute middles coordinates
P_121 = two_third * T.vertex(1,:) + one_third * T.vertex(2,:);
P_122 = one_third * T.vertex(1,:) + two_third * T.vertex(2,:);
P_131 = two_third * T.vertex(1,:) + one_third * T.vertex(3,:);
P_132 = one_third * T.vertex(1,:) + two_third * T.vertex(3,:);
P_231 = two_third * T.vertex(2,:) + one_third * T.vertex(3,:);
P_232 = one_third * T.vertex(2,:) + two_third * T.vertex(3,:);
P_123 = one_third * (T.vertex(1,:) + T.vertex(2,:) + T.vertex(3,:));

New_summit = compute_new_triangle_summit_coordinates(P_121, P_123, P_131);

V_new = [T.vertex(1,:); T.vertex(2,:); T.vertex(3,:); % already existing vertices
         P_121; P_122; P_131; P_132; P_231; P_232; P_123; New_summit]; % new vertices

% 11 = 9 - 1 + 3 new facets
% Sorted by top vertex first

F_new = [1 4 6;
         5 10 4;
         7 6 10;
         5 2 8;
         7 9 3;
         10 8 9;
         8 10 5;
         9 7 10;
         11 4 10;
         11 10 6;
         11 6 4];
 
% triangle / facet subdivision diagram
% New_summit is not visible here
% New_summit if the top of the tetrahedron which base is : [121 123 131]
%     
%                  1
%                /   \
%               /     \
%              /       \
%           121 _______ 131
%            /  \     /  \
%           /    \   /    \
%          /      \ /      \
%       122  ____ 123 _____ 132
%        / \      / \      / \
%       /   \    /   \    /   \
%      /     \  /     \  /     \
%     2 ____ 231 ____ 232 _____ 3

end


%% Split tetrahedron subfunction
function [V_new, F_new] = split_tetra(T)
%
% Input
%
% T : tetrahedron structure
%
% Outputs
%
% N : The six newly created vertices -coordinates-
% I : Four index line triplets of the four newly created facets
% P : The four previous vertices -coordinates-


% Compute middles coordinates
M_12 = 0.5 * (T.vertex(1,:) + T.vertex(2,:));
M_13 = 0.5 * (T.vertex(1,:) + T.vertex(3,:));
M_14 = 0.5 * (T.vertex(1,:) + T.vertex(4,:));
M_23 = 0.5 * (T.vertex(2,:) + T.vertex(3,:));
M_24 = 0.5 * (T.vertex(2,:) + T.vertex(4,:));
M_34 = 0.5 * (T.vertex(3,:) + T.vertex(4,:));

S_586 =  compute_new_triangle_summit_coordinates(M_12, M_23, M_13);
S_579 =  compute_new_triangle_summit_coordinates(M_12, M_14, M_24);
S_6107 = compute_new_triangle_summit_coordinates(M_13, M_34, M_14);
S_8910 = compute_new_triangle_summit_coordinates(M_23, M_24, M_34);

V_new = [T.vertex(1,:); T.vertex(2,:); T.vertex(3,:); T.vertex(4,:); % already existing vertices
         M_12; M_13; M_14; M_23; M_24; M_34; % new vertices
         S_586; S_579; S_6107; S_8910]; % new summits

% 12 = 4x3 new facets
% Sorted by trigo order per line

F_new = [11 5 8;
         11 8 6;
         11 6 5;...
         
         12 5 7;
         12 7 9;
         12 9 5;...
         
         13 6 10;
         13 10 7;
         13 7 6;...
         
         14 8 9;
         14 9 10;
         14 10 8;...
         
         1 5 6;
         1 6 7;
         1 7 5;...
         
         2 5 9;
         2 9 8;
         2 8 5;...
         
         3 6 8;
         3 8 10;
         3 10 6;...
                          
         4 7 10;
         4 10 9;
         4 9 7];
      
     
end


%% Compute new triangle summit coordinates subfunction
function [S] = compute_new_triangle_summit_coordinates(M, N, P)
%
% Input
%
% M = [Mx My Mz]
% N = [Nx Ny Nz]
% P = [Px Py Pz]
%
% in trigonometric / counter clockwise order
%
% Output
%
% S = [Sx Sy Sz] : new tetrahedron summit


a = norm(N-M);   % edge length
h = a*sqrt(2/3); % tetrahedron height

cross_prod = cross(P-N, M-N);
norm_cross_prod = cross_prod / norm(cross_prod);

S = mean([M; N; P], 1) - h * norm_cross_prod;

end


%% Triangles concatenation subfunction
function [V_array, T_array] = catriangles(Triangle_array)

S = size(Triangle_array,3);
V_array = zeros(3*S,3);
T_array = zeros(S,3);

for k = 1:S
    
    for i = 1:size(Triangle_array(:,:,k).vertex,1)
        
        V_array(3*(k-1)+i,:) = Triangle_array(:,:,k).vertex(i,:);
        
    end        
        
        a = Triangle_array(:,:,k).facet(1,1) + 3*(k-1);
        b = Triangle_array(:,:,k).facet(1,2) + 3*(k-1);
        c = Triangle_array(:,:,k).facet(1,3) + 3*(k-1);        
        
        T = [a b c];
        
        T_array(k,:) = T;
    
end

end


%% Display subfunction
function [] = disp_koch_snowflake(V, T, cmap)
%
% Input
%
% Triangle_array : Triangle array following dim 3

figure;
set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0]);
trisurf(T,V(:,1),V(:,2),V(:,3),'EdgeColor',cmap), shading interp, hold on;
colormap(cmap);
axis square, axis equal, axis tight, axis off;
grid off;
ax = gca;
ax.Clipping = 'off';
camlight left;
view(-90,19);

end


%% Remove duplicated vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)
tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);
end % remove_duplicated_vertices