function smoothedCurve = MAVG (curve,order)

curve = curve(:);
a = (order-1)/2;
curveEx = [flipud (curve(1:a));curve];

end