PRO REFRESH_IMAGE, Image=image, _Extra=extra, WID=wid, $
    yaxisSolarRadii=yaxisSolarRadii,imgTitle=imgTitle
  ;+
  ;
  ;   Name: REFRESH_IMAGE
  ;
  ;   Purpose: Refresh the image during the Contrast selection. Called by Xcolor
  ;
  ;-
    
  IF N_Elements(wid) NE 0 THEN WSet, wid
  plot_image, image, _Extra=extra, title=exist(imgTitle)?imgTitle:''
  IF exist(yaxisSolarRadii) THEN $
    axis,yaxis=1,yticks=5,ytickn=STRTRIM(string(float(yaxisSolarRadii), $
    FORMAT='(f5.2)'),1), ystyle=1,color=0,ytitle= 'Distance (Solar Radius)'
    
END