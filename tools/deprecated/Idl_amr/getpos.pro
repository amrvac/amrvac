FUNCTION GETPOS, ASPECT, POSITION=POSITION, MARGIN=MARGIN
;- Check arguments
if (n_params() ne 1) then $
  message, 'Usage: RESULT = GETPOS(ASPECT)'
if (n_elements(aspect) eq 0) then $
  message, 'Argument ASPECT is undefined'
;- Check keywords
if (n_elements(position) eq 0) then $
  position = [0.0, 0.0, 1.0, 1.0]
if (n_elements(margin) eq 0) then margin = 0.1
;- Get range limited aspect ratio and margin input values
aspect_val = (float(aspect[0]) > 0.01) < 100.0
margin_val = (float(margin[0]) > 0.0) < 0.49
;- Compute aspect ratio of position vector in this window
xsize = (position[2] - position[0]) * !d.x_vsize
ysize = (position[3] - position[1]) * !d.y_vsize
cur_aspect = ysize / xsize
;- Compute aspect ratio of this window
win_aspect = float(!d.y_vsize) / float(!d.x_vsize)
;- Compute height and width in normalized units
if (aspect_val ge cur_aspect) then begin
  height = (position[3] - position[1]) - 2.0 * margin
  width = height * (win_aspect / aspect_val)
endif else begin
  width = (position[2] - position[0]) - 2.0 * margin
  height = width * (aspect_val / win_aspect)
endelse
;- Compute and return position vector
xcenter = 0.5 * (position[0] + position[2])
ycenter = 0.5 * (position[1] + position[3])
x0 = xcenter - 0.5 * width
y0 = ycenter - 0.5 * height
x1 = xcenter + 0.5 * width
y1 = ycenter + 0.5 * height
return, [x0, y0, x1, y1]
END
