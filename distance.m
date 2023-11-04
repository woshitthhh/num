function d=distance(long1,lat1,long2,lat2)
    lat1=lat1*pi/180;
    lat2=lat2*pi/180;
    long1=long1*pi/180;
    long2=long2*pi/180;
    d=6378*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(long1-long2));
end
