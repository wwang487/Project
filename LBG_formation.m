function out_img = LBG_formation(orig_img)
 
orig_img = im2double(orig_img);

orig_height = size(orig_img,1);
    orig_width = size(orig_img,2);
    height = orig_height;
    width = orig_width;
    if 2*round(orig_height/2) ~= orig_height
        height = height + 1;
    end
    if 2*round(orig_width/2) ~= orig_width
        width = width + 1;
    end
    
    orig_img = imresize(orig_img,[height,width],'bilinear');
    gray_img = rgb2gray(orig_img);
    edge_img = edge(gray_img,'canny',0.25);
    
    LBG_img = rgbtoluv(orig_img);
    image_pixel = zeros(12,height*width/4);
        
    for i_node = 1 : height/2
        for j_node = 1 : width/2
            i = 2*i_node - 1;
            j = 2*j_node - 1;
            image_pixel(:,(i_node-1)*width/2 + j_node) = [LBG_img(i,j,1);LBG_img(i,j,2);LBG_img(i,j,3);...
                LBG_img(i+1,j,1);LBG_img(i+1,j,2);LBG_img(i+1,j,3);...
                LBG_img(i,j+1,1);LBG_img(i,j+1,2);LBG_img(i,j+1,3);...
                LBG_img(i+1,j+1,1);LBG_img(i+1,j+1,2);LBG_img(i+1,j+1,3)];
        end
    end
    
    LBG_pixel = LBG_VQ(image_pixel);
    LBG0_img = zeros(height,width,3);
    for j0 = 1 : height*width/4
        x_index = rem(j0-1, width/2)+1;
        y_index = round(2*(j0 - x_index)/width+1);
        x_index = 2 * x_index -1;
        y_index = 2 * y_index -1;
        
        LBG0_img(y_index,x_index,1)= LBG_pixel(1,j0);
        LBG0_img(y_index,x_index,2)= LBG_pixel(2,j0);
        LBG0_img(y_index,x_index,3)= LBG_pixel(3,j0);
        LBG0_img(y_index+1,x_index,1)= LBG_pixel(4,j0);
        LBG0_img(y_index+1,x_index,2)= LBG_pixel(5,j0);
        LBG0_img(y_index+1,x_index,3)= LBG_pixel(6,j0);
        LBG0_img(y_index,x_index+1,1)= LBG_pixel(7,j0);
        LBG0_img(y_index,x_index+1,2)= LBG_pixel(8,j0);
        LBG0_img(y_index,x_index+1,3)= LBG_pixel(9,j0);
        LBG0_img(y_index+1,x_index+1,1)= LBG_pixel(10,j0);
        LBG0_img(y_index+1,x_index+1,2)= LBG_pixel(11,j0);
        LBG0_img(y_index+1,x_index+1,3)= LBG_pixel(12,j0);
    end
    
 L_img = luvtorgb(LBG0_img);
 fR = L_img(:,:,1);
 fG = L_img(:,:,2);
 fB = L_img(:,:,3);
 sigma = 1;
 
 w = fspecial('average', 2);
 
 fR_filtered = imfilter(fR, w, 'replicate');      
 fG_filtered = imfilter(fG, w, 'replicate');      
 fB_filtered = imfilter(fB, w, 'replicate');
 fc_filtered = sigma*cat(3, fR_filtered+0.05, fG_filtered+0.1, fB_filtered);
 
 w1 = fspecial('laplacian',0.5);
 g1 = imfilter(fc_filtered,w1,'replicate');  
 fc_filtered =  fc_filtered -g1;  
 
 
 filtered_L_img = imadjust(fc_filtered,[0,1],[0.1,0.8]);
 
 %Gamma = 0.4;    
 %filtered_L_img = myExpEnhance(fc_filtered,Gamma); 
 filtered_L_img = remove_blur(filtered_L_img);

 for i = 1: height
     for j = 1:width
         if edge_img(i,j)~=0
             filtered_L_img(i,j,1)=0;
             filtered_L_img(i,j,2)=0;
             filtered_L_img(i,j,3)=0;
         end
     end
 end
    out_img = 0.8*filtered_L_img;
end
 

    
function LBG_Vector = LBG_VQ(input_pixel)
    LBG_Vector = input_pixel;
    pixel_num = size(input_pixel,2);
    
    training_times = 5;
    i = 1;
    Extra = [1;1;1;1;1;1;1;1;1;1;1;1];

    Lib_marker = ones(pixel_num,1);
    temp_marker = ones(pixel_num,1);
       
    while i <= training_times 
        marker_type = 2 ^ (i-1);        
        for c_val = 1 : marker_type            
            marker_indices = find(Lib_marker == c_val);
            C_sum = zeros(12,1);
            for j = 1:length(marker_indices)
                C_sum = C_sum + input_pixel(:,marker_indices(j,1));
            end
            C0 = round(C_sum./length(marker_indices));
            C1 = C0 - Extra;
            C2 = C0 + Extra;
            for k = 1:length(marker_indices)
                C = input_pixel(:,marker_indices(k,1));
                dC1 = sqrt(sum(C-C1).^2);
                dC2 = sqrt(sum(C-C2).^2);
                if dC1 <= dC2
                    temp_marker(marker_indices(k,1)) = 2 * c_val - 1;
                else 
                    temp_marker(marker_indices(k,1)) = 2 * c_val;
                end
            end            
        end
        Lib_marker = temp_marker;
        i = i+1;
    end
    
    for marker = 1:max(Lib_marker)
        Aver_index = find(Lib_marker == marker);
        Aver_size = size(Aver_index,1);
        Aver_matrix = zeros(12,Aver_size);
        for p = 1:Aver_size
            Aver_matrix(:,p) = input_pixel(:,Aver_index(p,1));
        end
        
        C0 = median(Aver_matrix,2); 
        C_Lightness = zeros(4,5);
        for i = 1:4
            C_min = min(Aver_matrix(3*i -2,:));
            C_max = max(Aver_matrix(3*i -2,:));
            for j = 1:5
                C_Lightness(i,j) = round(C_max - (5-j)*(C_max - C_min)/4);
            end
        end
        
        for q = 1:length(Aver_index)
            for i = 1:4
                min_dist = LBG_Vector(3*i-2, Aver_index(q))-C_Lightness(i,1);
                index = 1;
                for j = 2:5
                    if LBG_Vector(3*i-2, Aver_index(q))-C_Lightness(i,j) < min_dist
                        index = j;
                        min_dist = LBG_Vector(3*i-2, Aver_index(q))-C_Lightness(i,j);
                    end
                end
                LBG_Vector(3*i-2, Aver_index(q))= C_Lightness(i,index)-3;
            end
        end        
    end          
end



function dst_img=myExpEnhance(src_img,Gamma)    
    src_img = mat2gray(src_img,[0 255]);
    C = 1;    
    g2 = C*(src_img.^Gamma);   
    max=255;  
    min=0;  
    dst_img=uint8(g2*(max-min)+min); 
end


