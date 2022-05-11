function drawVessels(Vessel)
%drawVesselsNew - Draws a vascular tree using line segments. Colours segments to correspond to vasculature catergory. (NOTE: Not reccomended if vascular tree contains >20,000 nodes)
%
% Syntax:  drawVesselsNew(Vessel)
%
% Inputs:
%    Vessel - The vessel tree for drawing. Required to be a Nx7 matrix
%    where N is the number of nodes on the vessel tree. Spatial coordinates
%    are in columns 3, 4, & 5 for [x y z] coordinate and column 7 contains
%    the connection node for the vessel tree.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Luke Fulford - luke.fulford@ed.ac.uk
% Date Modified: 11/01/2019

grid on % Turns grid on.
hold on
view([-45 30])


for Pnt=2:size(Vessel,1) % For all line segments in vessel tree.
    
    PntX=Vessel(Pnt,3); PntY=Vessel(Pnt,4); PntZ=Vessel(Pnt,5);
    % Establish current point locations.
    
    PntConn=Vessel(Pnt,7);
    PntConnX=Vessel(PntConn,3); PntConnY=Vessel(PntConn,4); PntConnZ=Vessel(PntConn,5);
    % Establish connecting point locations.
    
    %if (Vessel(Pnt,2)==2) || (Vessel(Pnt,2)==12) % 2 is LACA
    if mod(Vessel(Pnt,2),10)==2
        ColourString='g';
        %plot3([PntConnY PntY],[PntConnX PntX],[PntConnZ PntZ],'Color',ColourString);
    elseif mod(Vessel(Pnt,2),10)==3 %(Vessel(Pnt,2)==3) || (Vessel(Pnt,2)==13) % 3 is LMCA
        ColourString='c';
        %plot3([PntConnY PntY],[PntConnX PntX],[PntConnZ PntZ],'Color',ColourString);
    elseif mod(Vessel(Pnt,2),10)==4 %(Vessel(Pnt,2)==4) || (Vessel(Pnt,2)==14) % 4 is RMCA
        ColourString='r';
    elseif mod(Vessel(Pnt,2),10)==5 %(Vessel(Pnt,2)==5) || (Vessel(Pnt,2)==15) % 5 is RACA
        ColourString='y';
    elseif mod(Vessel(Pnt,2),10)==6 %(Vessel(Pnt,2)==6) || (Vessel(Pnt,2)==16) % 6 is LPCA
        ColourString='m';
    elseif mod(Vessel(Pnt,2),10)==7 %(Vessel(Pnt,2)==7) || (Vessel(Pnt,2)==17) % 7 is RPCA
        ColourString='b';
    elseif mod(Vessel(Pnt,2),10)==0 %Vessel(Pnt,2)==10
        ColourString=[0.9290, 0.6940, 0.1250];
    else
        ColourString='k';
    end
    % This colours the vessels according to vessel classification. This
    % defaults to black for vessels from the orginal vessel tree and red
    % for vessels generated or split from vessel generation.
    
    plot3([PntConnY PntY],[PntConnX PntX],[PntConnZ PntZ],'Color',ColourString, 'LineWidth', 2);
    % Plot vessel lines.
    
    drawnow limitrate nocallbacks % Updates figure as it loops.
end
drawnow % Finalises figure.
