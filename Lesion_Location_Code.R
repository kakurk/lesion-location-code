F_ROI_lesion_stat3 <- function(IP_subject,IP_Lesion_file_l,IP_ROI_coordinate,IP_output_file_tag,IP_ROI_IDs,IP_case_names)#the main funtion begins
{
  
  #Load the fmri library
  library("fmri")
  
  #----------------------------------------------------------------------------------------------------------------
  
  #which positions on the IP_Lesion_file_l belong to the IP subject
  indx <- which(IP_case_names == IP_subject)
  
  #----------------------------------------------------------------------------------------------------------------
  
  #read the lesion file(s) and make an overlap mask
  print(indx)
  print(IP_Lesion_file_l[indx])
  
  #read the first file
  lesion_mask <- read.NIFTI(IP_Lesion_file_l[indx[1]], level = 0.75,setmask=FALSE)
  overlap_mask_data <- extract.data(lesion_mask, what = "data")
  
  #now overlap the other masks on the first mask
  c1 <- 1
  while(c1 < length(indx))
  {
    lesion_mask <- read.NIFTI(IP_Lesion_file_l[indx[c1]], level = 0.75, setmask=FALSE)
    lesion_mask_data <- extract.data(lesion_mask, what = "data")
    
    overlap_mask_data <- overlap_mask_data + lesion_mask_data
    
    c1 <- c1 + 1
  }
  
  # indices of all voxels > 0
  voxels <- which(overlap_mask_data > 0, arr.ind = TRUE)
  
  # remove the fourth column
  voxels <- voxels[,1:3]
  
  #write.ANALYZE(array(overlap_mask_data,c(lesion_mask$dim[1],lesion_mask$dim[2],lesion_mask$dim[3],1)),file=paste("anatomic",sep=""))
  #lesion_threshold <- 0.1*max(overlap_mask_data) #please note due to normalization of lesion mask the lesion map is smudged a bit
  #print(c(max(overlap_mask_data),lesion_threshold))
  
  #-----------------------------------------------------------------------------------------------------------
  
  #count total lesion voxels
  lesion_counter <- nrow(voxels)
  
  #write!
  write(c(IP_subject,lesion_counter),file=paste("Total_lesion_voxels",IP_subject,".dat",sep=""),ncol=2000,append=TRUE,sep =" ")
  
  #-----------------------------------------------------------------------------------------------------------
  
  #first read powers coordinate file and then evaluate the overlap between lesion and the ROIs
  p_cood <- read.table(IP_ROI_coordinate, header=TRUE)
  
  # initalize a vector to collect distances
  distance_distribution <- vector(mode = "numeric", length = nrow(voxels))
  
  #now iterate through all given ROI_IDs and find out the mean distance from lesion voxel
  for(c2 in 1:length(IP_ROI_IDs))#for.c2.begins  #PARA
  {
    
    #fetch the roi_id you want to analyze
    roi_id <- IP_ROI_IDs[c2]
    
    #find the voxel closest to this power's coordinate
    #if you set the mm boxes to 0, 0, 0, you will get voxel space as 91, 107, 89
    vox_x_cood <- round((91 - p_cood$x[roi_id]))
    vox_y_cood <- round((107 + p_cood$y[roi_id])) 
    vox_z_cood <- round((89 + p_cood$z[roi_id])) 
    the_coordinates_of_this_power264_voxel <- c(vox_x_cood, vox_y_cood, vox_z_cood)
    
    # lets loop over the voxels we found
    for(v in 1:nrow(voxels)){
      
      # the coordinates (x,y,z) of the current voxel found in the lesion mask
      the_coordinates_of_this_voxel <- as.numeric(voxels[v, 1:3])
      
      # the distance between this lesion mask voxel and this power 264 coordinate 
      temp_distance <- dist(rbind(the_coordinates_of_this_power264_voxel, the_coordinates_of_this_voxel), method="euclidean")
      
      # concatenate into a list
      distance_distribution[v] <- temp_distance

    }
    
    # print mean, median, and minimum distances to the console
    print(c(c2, mean(distance_distribution), median(distance_distribution), min(distance_distribution)))
    
    # write mean, median, and minimum distances to a file
    write(c(c2, mean(distance_distribution), median(distance_distribution), min(distance_distribution)), file=IP_output_file_tag, ncol=2000, append=TRUE, sep =" ")
    
  }#for.c2.ends
  
}#the main funtion ends

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


#lesion_masks

l_sub1a <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.020/nnew_cot2fl3dtrap2swis007a1001-label.nii"
l_sub1b <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.020/nnew_T2FLAIRAXIALs005a1001-label.nii"

l_sub2a <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.023/nnew_20110127_104146t2fl3dtrap2swis011a1001-label.nii"
l_sub2b <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.023/nnew_20110127_104146T2FLAIRAXIALs009a1001-label.nii"

l_sub3a <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.026/nnew_20111201_095347t2fl3dtrap2swis017a1001-label.nii"
l_sub3b <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.026/nnew_20111201_095347T2FLAIRAXIALs015a1001-label.nii"

l_sub4a <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.027/nnew_20120222_180702t2fl3dtrap2swis019a1001-label.nii"

l_sub5a <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.030/nnew_20120620_171657T1MPRAGEIsos002a1001-label.nii"
l_sub5b <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.030/nnew_20120620_171657t2fl3dtrap2swis017a1001-label.nii"
l_sub5c <- "/Users/kylekurkela/Desktop/HillaryLab/long.1.030/nnew_20120620_171657T2FLAIRAXIALs015a1001-label.nii"


Lesion_file_l <-  c(l_sub1a, l_sub1b, l_sub2a, l_sub2b, l_sub3a, l_sub3b, l_sub4a, l_sub5a, l_sub5b, l_sub5c)

case_names <- c(1, 1, 2, 2, 3, 3, 4, 5, 5, 5)

#---------------------------------------------------------------------------------------------
#remember to change the ROI_coordinate path to account for the MNI coordinate information

if(length(case_names) == length(Lesion_file_l))
{
  for(c1 in 1:5)
  {
    temp_indx <- which(case_names == c1) #this will return you a list
    print(temp_indx)
    if(Lesion_file_l[temp_indx[1]] != "") #access the first position of this list and see if it is "" or not
    {
      ROI_coordinate <- '/Users/kylekurkela/Desktop/HillaryLab/powers_MNI_coordinate.dat' 
      Output_file_tag <- paste('./Report_case_',c1,sep="") 
      ROI_IDs <- c(1:264)
      F_ROI_lesion_stat3(c1,Lesion_file_l,ROI_coordinate,Output_file_tag,ROI_IDs,case_names)
    }
  }
}

if(length(case_names) != length(Lesion_file_l))
{
  print("length(case_names) is not equal to length(Lesion_file_l).") 
}