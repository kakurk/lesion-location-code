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
    lesion_mask <- read.NIFTI(IP_Lesion_file_l[indx[c1]], level = 0.75,setmask=FALSE)
    lesion_mask_data <- extract.data(lesion_mask, what = "data")
    
    overlap_mask_data <- overlap_mask_data + lesion_mask_data
    
    c1 <- c1 + 1
  }
  
  
  #write.ANALYZE(array(overlap_mask_data,c(lesion_mask$dim[1],lesion_mask$dim[2],lesion_mask$dim[3],1)),file=paste("anatomic",sep=""))
  #lesion_threshold <- 0.1*max(overlap_mask_data) #please note due to normalization of lesion mask the lesion map is smudged a bit
  #print(c(max(overlap_mask_data),lesion_threshold))
  
  #-----------------------------------------------------------------------------------------------------------
  
  #count total lesion voxels
  
  
  lesion_counter <- 0
  
  for(c3x in 1:lesion_mask$dim[1])#for.c3x.begins
  {
    for(c3y in 1:lesion_mask$dim[2])#for.c3y.begins
    {
      for(c3z in 1:lesion_mask$dim[3])#for.c3z.begins
      {
        
        lesion_val <- overlap_mask_data[c3x,c3y,c3z,1]
        if(lesion_val > 0)
        {lesion_counter <- lesion_counter + 1}
        
        
      }#for.c3z.ends
    }#for.c3y.ends
  }#for.c3x.ends
  
  
  write(c(IP_subject,lesion_counter),file=paste("Total_lesion_voxels",IP_subject,".dat",sep=""),ncol=2000,append=TRUE,sep =" ")
  
  #-----------------------------------------------------------------------------------------------------------
  
  #first read powers coordinate file and then evaluate the overlap between lesion and the ROIs
  
  p_cood <- read.table(IP_ROI_coordinate,header=TRUE)
  
  
  
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
    
    distance_distribution <- c()
    
    #now calculate the distance of ROI center to all lesion voxels
    for(c3x in 1:lesion_mask$dim[1])#for.c3x.begins
    {
      for(c3y in 1:lesion_mask$dim[2])#for.c3y.begins
      {
        for(c3z in 1:lesion_mask$dim[3])#for.c3z.begins
        {
          
          lesion_val <- overlap_mask_data[c3x,c3y,c3z,1]
          if(lesion_val > 0)
          { 
            temp_distance <- dist(rbind(c(vox_x_cood,vox_y_cood,vox_z_cood),c(c3x,c3y,c3z)),method="euclidean")
            distance_distribution <- c(distance_distribution,temp_distance)
          }
          
          
        }#for.c3z.ends
      }#for.c3y.ends
    }#for.c3x.ends
    
    print(c(c2,mean(distance_distribution),median(distance_distribution), min(distance_distribution)))
    
    write(c(c2,mean(distance_distribution),median(distance_distribution), min(distance_distribution)),file=IP_output_file_tag,ncol=2000,append=TRUE,sep =" ")
    
  }#for.c2.ends
  
  
  #----------------------------------------------------------------------------------------------------------------
  
}#the main funtion ends

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


#lesion_masks


l_sub1a <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.020/nnew_cot2fl3dtrap2swis007a1001-label.nii"
l_sub1b <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.020/nnew_T2FLAIRAXIALs005a1001-label.nii"

l_sub2a <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.023/nnew_20110127_104146t2fl3dtrap2swis011a1001-label.nii"
l_sub2b <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.023/nnew_20110127_104146T2FLAIRAXIALs009a1001-label.nii"

l_sub3a <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.026/nnew_20111201_095347t2fl3dtrap2swis017a1001-label.nii"
l_sub3b <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.026/nnew_20111201_095347T2FLAIRAXIALs015a1001-label.nii"

l_sub4a <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.027/nnew_20120222_180702t2fl3dtrap2swis019a1001-label.nii"

l_sub5a <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.030/nnew_20120620_171657T1MPRAGEIsos002a1001-label.nii"
l_sub5b <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.030/nnew_20120620_171657t2fl3dtrap2swis017a1001-label.nii"
l_sub5c <- "/storage/home/epg5130/work/lesion_distance_analysis/long.1.030/nnew_20120620_171657T2FLAIRAXIALs015a1001-label.nii"


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
      ROI_coordinate <- '/gpfs/group/nad12/legacy/fgh3/Studies/Lesions/Lesion_library/TBI Complete/powers_MNI_coordinate.dat' 
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
