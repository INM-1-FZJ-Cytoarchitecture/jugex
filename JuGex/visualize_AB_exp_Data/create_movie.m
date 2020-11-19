patch_hdl_fp1.Visible='Off';
patch_hdl_fp2.Visible='Off';

area1_name='FP2';
area2_name='FP1';

h_bubbles=findobj(h_curfig,'Type','surface','-not','Tag','zslice','-not','Tag','xslice','-not','Tag','yslice');
set(h_bubbles,'Visible','Off')
h_crosses=findobj(h_curfig,'Type','patch','-not','Tag','brain_surface');
set(h_crosses,'Visible','Off')

b=repmat(shift_between_map_and_ref,size(area1_coords,2),1);
area1_coords=area1_coords-int32(b)';
mean_coords1=mean(area1_coords,2);
b=repmat(shift_between_map_and_ref,size(area2_coords,2),1);
area2_coords=area2_coords-int32(b)';
mean_coords2=mean(area2_coords,2);

h_arrow_area1=arrow3([190,120,80],[190,100,80],'r');
h_text_area1=text(190,120,100,area1_name,'Color','red','FontSize',14);

h_arrow_area2=arrow3([200,40,80],[200,60,80],'b');
h_text_area2=text(200,40,80,area2_name,'Color','blue','FontSize',14);


F(1) = getframe(h_figure);
frame_counter=2;
capture_movie=1;

if capture_movie==1
    patch_hdl_fp1.Visible='On';
patch_hdl_fp2.Visible='On';
    for i=0:0.01:0.3
patch_hdl_fp1.FaceAlpha	=i;
patch_hdl_fp2.FaceAlpha	=i;
[F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
    for i=0.3:-0.01:0
patch_hdl_fp1.FaceAlpha	=i;
patch_hdl_fp2.FaceAlpha	=i;
[F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
    patch_hdl_fp1.Visible='Off';
patch_hdl_fp2.Visible='Off';

end

if capture_movie==1
h_bubbles=findobj(h_curfig,'Type','surface','-not','Tag','zslice','-not','Tag','xslice','-not','Tag','yslice');
set(h_bubbles,'Visible','On')
h_text_area1.Visible='Off';
h_text2_area1=text(190,192,62,{'Spheres represent ',['tissue blocks in ' area1_name '.' ],'Colors indicate', 'expression levels' '(z-scores).'},'Color','red','FontSize',12,'BackgroundColor','w','EdgeColor','r');
[F,frame_counter] = capture_n_frames( 120,frame_counter,F,h_figure );
h_crosses=findobj(h_curfig,'Type','patch','-not','Tag','brain_surface');
set(h_crosses,'Visible','On')
h_text_area2.Visible='Off';
h_text3_area1=text(200,40,80,{['Crosses represent '],[' tissue blocks in ' area2_name '.' ]},'Color','blue','FontSize',12,'BackgroundColor','w','EdgeColor','b');
[F,frame_counter] = capture_n_frames( 60,frame_counter,F,h_figure );

 delete(h_text2_area1);
 delete(h_text3_area1);
end

 h_text_area2.Visible='Off';
 h_text_area1.Visible='Off';
  set( h_arrow_area1,'Visible','Off')
   set( h_arrow_area2,'Visible','Off')
%clear h_arrow_area1 h_arrow_area2 h_text_area1 h_text_area1
for zoom_f=0:100
    frame_counter;
%     zzz=1+(zoom_f/1000);
%     zoom(zzz)
    zzz=1+(zoom_f/10000);
    camzoom(zzz);
    drawnow
    if capture_movie==1
        [F,frame_counter] = capture_n_frames( 1,frame_counter,F,h_figure);
    end
end
%%%% try movie
%frame_counter=2;
[az,el] = view;

for angel=az:2:-50
    frame_counter;
    view(angel,14);
    %lightsource=camlight(lightsource,10,180);
    lightsource=camlight(lightsource,'headlight'); 
    drawnow
    if capture_movie==1
        [F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
end



for angel=-50:-2:-140
    frame_counter;
    view(angel,14);
    %lightsource=camlight(lightsource,10,180);
    lightsource=camlight(lightsource,'headlight'); 
    drawnow
    if capture_movie==1
        [F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
end

top_angel=[14:2:80];
top_angel=[top_angel,fliplr(top_angel)];
top_angel(69:91)=14;
tac=1;
for angel=-140:-50
    frame_counter;
    view(angel,top_angel(tac));
    %set(gca,'xlim',[tac 250]);
    if angel<80||angel>100
        %lightsource=camlight(lightsource,10,180);
        lightsource=camlight(lightsource,'headlight'); 
    end
    drawnow
    if capture_movie==1
    [F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
    tac=tac+1;
end

% hier könnten die karten eingeblendet werden
% schritte davor kürzer
patch_hdl_fp1.Visible='On';
patch_hdl_fp2.Visible='On';

for zoom_f=0:10
    frame_counter;
    zzz=1+(zoom_f/100);
    camzoom(zzz);
    drawnow
    if capture_movie==1
    [F,frame_counter] = capture_n_frames( 4,frame_counter,F,h_figure );
    end
end

if capture_movie==1
    for i=0:0.01:0.3
patch_hdl_fp1.FaceAlpha	=i;
patch_hdl_fp2.FaceAlpha	=i;
[F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
    for i=0.3:-0.01:0.1
patch_hdl_fp1.FaceAlpha	=i;
patch_hdl_fp2.FaceAlpha	=i;
[F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
end

%2ter zoom




patch_hdl_fp1.Visible='Off';
patch_hdl_fp2.Visible='Off';

    
for zoom_f=0:10
    frame_counter;
    zzz=1-(zoom_f/100);
    camzoom(zzz)
    drawnow
    if capture_movie==1
        [F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
end    
 
for zoom_f=0:10
    frame_counter;
    zzz=1-(zoom_f/100);
    camzoom(zzz)
    drawnow
    if capture_movie==1
        [F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
end 



stepis=linspace(-50,az,30);
for angel=1:30
    frame_counter;
    view(stepis(angel),14);
   % lightsource=camlight(lightsource,10,180);
   lightsource=camlight(lightsource,'headlight'); 
    drawnow
    if capture_movie==1
        [F,frame_counter] = capture_n_frames( 2,frame_counter,F,h_figure );
    end
end


    
disp('ready. Still have to save!');
name_avi='test.avi';
writerObj = VideoWriter(name_avi);
writerObj.FrameRate = 30;
writerObj.Quality = 100;
open(writerObj);
writeVideo(writerObj,F);
close(writerObj);
disp([name_avi ' saved!']);