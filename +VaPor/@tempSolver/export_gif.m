function export_gif(obj, filename, clims)
h = figure;
set(gca, 'FontName', 'Fira Sans Light')
axis tight manual % this ensures that getframe() returns a consistent size
%filename = 'variableCore.gif';

for n = 1:size(obj.Transient_Temp.tissue ,2)
    fprintf('Processing Frame: %d \n', n);
    % Marshall the data from the T_TransientStore into a format we can use
    % to plot the relevant image:
    %%%%%% Reshaping Results from Solver %%%%%%%%%%%%
        obj.plot_trans_brain(n, 30);
        drawnow;
        % Clamp the caxis:
        c = colorbar;
        c.FontName = 'Fira Sans Light';
        caxis(clims);
        colormap jet
        shading interp;
        frameTitle = sprintf('Brain temperature after %.1f mins', (n-1)*obj.save_timestep/60);
        title(frameTitle, 'FontName', 'Fira Sans Light');
        axis off;
    
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.01);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.01);
    end
end
end