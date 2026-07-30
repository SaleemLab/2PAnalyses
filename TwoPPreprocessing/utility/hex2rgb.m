function rgb = hex2rgb(hexStr)
    hexStr = strrep(hexStr, '#', '');
    rgb = [hex2dec(hexStr(1:2)), hex2dec(hexStr(3:4)), hex2dec(hexStr(5:6))] / 255;
end