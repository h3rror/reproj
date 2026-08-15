function make_snapshot_movie(Xsnap, Ny, Nx, filename, fps,skip)
% MAKE_SNAPSHOT_MOVIE  Render vectorized snapshots as a movie.
%
%   make_snapshot_movie(Xsnap, Ny, Nx, filename, fps)
%
%   Xsnap    : Ndof x Nsnap matrix, each column a vectorized 2D field
%              (Ndof = Ny*Nx, consistent with C(:) / reshape(.,Ny,Nx))
%   Ny, Nx   : grid dimensions used to reshape each column
%   filename : output video file, e.g. 'advection.mp4'
%   fps      : frames per second (default 15)
%
% Example:
%   make_snapshot_movie(Xsnap, Ny, Nx, 'advection.mp4', 20);

if nargin < 5, fps = 15; end

v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = fps;
open(v);

clim = [min(Xsnap(:)), max(Xsnap(:))];   % fixed color scale across frames

fig = figure('Color','w');
for k = 1:skip:size(Xsnap,2)
    frame_field = reshape(Xsnap(:,k), Ny, Nx);
    imagesc(frame_field, clim);
    axis equal tight; colorbar;
    title(sprintf('Snapshot %d / %d', k, size(Xsnap,2)));
    drawnow;
    writeVideo(v, getframe(fig));
end

close(v);
close(fig);
fprintf('Movie saved to %s (%d frames @ %d fps)\n', filename, size(Xsnap,2), fps);
end