function roll_mean(img; voxel_x = 1, voxel_y = 1, voxel_z = 100, boundary_mode ="replicate")
     if voxel_x == 0
          voxel_x = size(img,1)
     end

     if voxel_y == 0
          voxel_y == size(img,2)
     end

     kernel = centered(ones(voxel_x,voxel_y,voxel_z) / (voxel_x*voxel_y*voxel_z))
     imfilter(img, kernel, boundary_mode)
end

function deltaF(img; mode = :mean, kwargs...)
     background = roll_mean(img; kwargs...)
     img .- background
end

function deltaF_F(img; mode = :mean, kwargs...)
     if mode == :rolling_mean
          background = roll_mean(img; kwargs...)
     elseif mode == :mean
          background = mean(img, dims = 3)[:,:,1]
     end
     (img .- background)./background
end