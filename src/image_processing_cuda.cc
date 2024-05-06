

// __global__ void computeTangentKernel(cv::cuda::PtrStepSz<cv::Vec3f> src_grad, cv::cuda::PtrStepSz<cv::Vec3f> dst_tangent, int w, int h) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i >= w || j >= h)
//         return;

//     cv::Vec3f src_grad_pixel = src_grad(j, i);
//     cv::Vec3f dst_tangent_pixel;

//     if (src_grad_pixel[2] == 0.0f) {
//         // Assuming vzero sets vector to zero
//         dst_tangent_pixel = cv::Vec3f(0.0f, 0.0f, 0.0f);
//     } else {
//         // Compute gradient vector ignoring the z component
//         float gx = src_grad_pixel[0];
//         float gy = src_grad_pixel[1];

//         // Axis z is predefined as (0, 0, 1)
//         float ax = 0.0f;
//         float ay = 0.0f;
//         float az = 1.0f;

//         // Compute cross product (gradient x axis_z)
//         dst_tangent_pixel[0] = gy * az - 0.0f; // gy * 1 - 0
//         dst_tangent_pixel[1] = 0.0f - gx * az; // 0 - gx * 1
//         dst_tangent_pixel[2] = gx * ay - gy * ax; // gx * 0 - gy * 0
//     }

//     dst_tangent(j, i) = dst_tangent_pixel;
// }

// void get_tangent(cv::cuda::GpuMat& src_grad, cv::cuda::GpuMat& dst_tangent, int w, int h) {
//     dim3 threadsPerBlock(16, 16);
//     dim3 numBlocks((w + threadsPerBlock.x - 1) / threadsPerBlock.x,
//                    (h + threadsPerBlock.y - 1) / threadsPerBlock.y);

//     computeTangentKernel<<<numBlocks, threadsPerBlock>>>(src_grad, dst_tangent, w, h);
// }
