package gamemobile.kmu.ac.kr.vrapp3;

import android.content.Context;
import android.content.res.AssetManager;
import android.graphics.Bitmap;
import android.opengl.GLES10;
import android.opengl.GLES20;
import android.opengl.GLUtils;
import android.opengl.Matrix;
import android.os.Environment;
import android.util.Log;

import com.google.vr.sdk.base.Eye;
import com.google.vr.sdk.base.GvrView;
import com.google.vr.sdk.base.HeadTransform;
import com.google.vr.sdk.base.Viewport;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;

import javax.microedition.khronos.egl.EGLConfig;
import javax.microedition.khronos.opengles.GL10;

/**
 * Created by user on 2018-03-22.
 */

public class VRApp3Renderer implements GvrView.StereoRenderer { // GvrView.Renderer {
    private int mPositionHandle;
    private int mTexCoordinateHandle;
    private int mTextureUniformHandle;
    private int mTextures[];

    private int mProgramHandle;

    private final int mBytesPerFloat = 4;
    /** Size of the position data in elements. */
    private final int mPositionDataSize = 3;
    private final int mTexCoordinateDataSize = 2;

    /** Store our model data in a float buffer. */
    private final FloatBuffer mScreenPosition;
    private final FloatBuffer mTextureCoordinate;

    int texW = 640;//canvas.getWidth();
    int texH = 480;//canvas.getHeight();

    private int inv = 0;
    private final boolean bFile = false;

    /** The buffer holding the vertices */
    private FloatBuffer vertexBuffer;
    /** The buffer holding the texture coordinates */
    private FloatBuffer textureBuffer;

    Context context;
    float camOrigTarg[];

    final float[] screenPosition =
        {
                -1.0f, 1.0f, 0.0f,
                -1.0f,-1.0f, 0.0f,
                1.0f, 1.0f, 0.0f,
                -1.0f,-1.0f, 0.0f,
                1.0f,-1.0f, 0.0f,
                1.0f, 1.0f, 0.0f
        };


    final float[] textureCoordinateData =
        {
                // Front face
                0.0f, 0.0f,
                0.0f, 1.0f,
                1.0f, 0.0f,
                0.0f, 1.0f,
                1.0f, 1.0f,
                1.0f, 0.0f
        };

    public native String stringFromJNI();
    public native float[] initSmallPtGPU(int u, int f, String k, int w, int h, String s, String r, AssetManager asset);
    public native int[] updateRendering();
    public native void finishRendering();
    public native void touchFunc(int deltax, int deltay);
    //Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_touchFunc
    public native void reinitCamera(float origx, float origy, float origz, float targx, float targy, float targz);

    // Used to load the 'native-lib' library on application startup.
    static
    {
        System.loadLibrary("native-lib");
    }

    public VRApp3Renderer(Context context) {
        this.context = context;

        // Initialize the buffers.
        mScreenPosition = ByteBuffer.allocateDirect(screenPosition.length * mBytesPerFloat).order(ByteOrder.nativeOrder()).asFloatBuffer();
        mTextureCoordinate = ByteBuffer.allocateDirect(textureCoordinateData.length * mBytesPerFloat).order(ByteOrder.nativeOrder()).asFloatBuffer();

        mScreenPosition.put(screenPosition).position(0);
        mTextureCoordinate.put(textureCoordinateData).position(0);
    }

    private void saveToFile(Bitmap leftBmp, Bitmap rightBmp) {
        String strFN;
        FileOutputStream out = null;

        if (leftBmp != null) {
            strFN = Environment.getExternalStorageDirectory() + "/" + context.getPackageName() + "/images/left_image_" + new Integer(inv).toString() + ".png";

            try {
                out = new FileOutputStream(strFN);
                leftBmp.compress(Bitmap.CompressFormat.PNG, 100, out); // bmp is your Bitmap instance
                // PNG is a lossless format, the compression factor (100) is ignored
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                try {
                    if (out != null) {
                        out.close();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        if (rightBmp!=null) {
            strFN = Environment.getExternalStorageDirectory() + "/" + context.getPackageName() + "/images/right_image_" + new Integer(inv).toString() + ".png";

            try {
                out = new FileOutputStream(strFN);
                rightBmp.compress(Bitmap.CompressFormat.PNG, 100, out); // bmp is your Bitmap instance
                // PNG is a lossless format, the compression factor (100) is ignored
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                try {
                    if (out != null) {
                        out.close();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    @Override
    public void onSurfaceCreated(EGLConfig eglConfig) {
        Log.e("VRApp3Renderer", "Start of onSurfaceCreated");

        final String vertexShader =
                "attribute vec4 a_Position;            \n"		// Per-vertex position information we will pass in.
                + "attribute vec2 a_TexCoordinate;       \n"		// Per-vertex position information we will pass in.
                + "varying vec2 v_TexCoordinate;         \n"
                + "                                      \n"
                + "void main()                           \n"		// The entry point for our vertex shader.
                + "{                                     \n"
                + "   v_TexCoordinate = a_TexCoordinate; \n"
                + "   gl_Position = a_Position;          \n" 	// gl_Position is a special variable used to store the final position.
                + "                                  	   \n"     // Multiply the vertex by the matrix to get the final point in
                + "}                                     \n";    // normalized screen coordinates.

        final String fragmentShader =
                "precision mediump float;					\n"
                + "                                      \n"
                + "uniform sampler2D u_Texture;          \n" // The input texture.
                + "varying vec2 v_TexCoordinate;         \n"
                + "              								\n"
                + "void main()              				\n"
                + "  {											\n"
                + "     gl_FragColor = texture2D(u_Texture, v_TexCoordinate);\n"
                + "  }              							\n";

        // Load in the vertex shader.
        int vertexShaderHandle = GLES20.glCreateShader(GLES20.GL_VERTEX_SHADER);
        if (vertexShaderHandle != 0)
        {
            // Pass in the shader source.
            GLES20.glShaderSource(vertexShaderHandle, vertexShader);

            // Compile the shader.
            GLES20.glCompileShader(vertexShaderHandle);

            // Get the compilation status.
            final int[] compileStatus = new int[1];
            GLES20.glGetShaderiv(vertexShaderHandle, GLES20.GL_COMPILE_STATUS, compileStatus, 0);

            // If the compilation failed, delete the shader.
            if (compileStatus[0] == 0)
            {
                GLES20.glDeleteShader(vertexShaderHandle);
                vertexShaderHandle = 0;
            }
        }
        if (vertexShaderHandle == 0)
        {
            throw new RuntimeException("Error creating vertex shader.");
        }

        // Load in the fragment shader shader.
        int fragmentShaderHandle = GLES20.glCreateShader(GLES20.GL_FRAGMENT_SHADER);
        if (fragmentShaderHandle != 0)
        {
            // Pass in the shader source.
            GLES20.glShaderSource(fragmentShaderHandle, fragmentShader);

            // Compile the shader.
            GLES20.glCompileShader(fragmentShaderHandle);

            // Get the compilation status.
            final int[] compileStatus = new int[1];
            GLES20.glGetShaderiv(fragmentShaderHandle, GLES20.GL_COMPILE_STATUS, compileStatus, 0);

            // If the compilation failed, delete the shader.
            if (compileStatus[0] == 0)
            {
                GLES20.glDeleteShader(fragmentShaderHandle);
                fragmentShaderHandle = 0;
            }
        }
        if (fragmentShaderHandle == 0)
        {
            throw new RuntimeException("Error creating fragment shader.");
        }

        // Create a program object and store the handle to it.
        mProgramHandle = GLES20.glCreateProgram();
        if (mProgramHandle != 0)
        {
            // Bind the vertex shader to the program.
            GLES20.glAttachShader(mProgramHandle, vertexShaderHandle);

            // Bind the fragment shader to the program.
            GLES20.glAttachShader(mProgramHandle, fragmentShaderHandle);

            // Bind attributes
            GLES20.glBindAttribLocation(mProgramHandle, 0, "a_Position");

            // Link the two shaders together into a program.
            GLES20.glLinkProgram(mProgramHandle);

            // Get the link status.
            final int[] linkStatus = new int[1];
            GLES20.glGetProgramiv(mProgramHandle, GLES20.GL_LINK_STATUS, linkStatus, 0);

            // If the link failed, delete the program.
            if (linkStatus[0] == 0)
            {
                GLES20.glDeleteProgram(mProgramHandle);
                mProgramHandle = 0;
            }
        }
        if (mProgramHandle == 0)
        {
            throw new RuntimeException("Error creating program.");
        }

        mPositionHandle = GLES20.glGetAttribLocation(mProgramHandle, "a_Position");
        mTexCoordinateHandle = GLES20.glGetAttribLocation(mProgramHandle, "a_TexCoordinate");
        mTextureUniformHandle = GLES20.glGetUniformLocation(mProgramHandle, "u_Texture");

        mTextures = new int[2];

        GLES20.glGenTextures(2, mTextures, 0);
        if (mTextures[0] == 0 || mTextures[1] == 0)
        {
            throw new RuntimeException("Error loading texture.");
        }

        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, mTextures[0]);

        GLES20.glTexParameterf(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MIN_FILTER, GLES20.GL_LINEAR);
        GLES20.glTexParameterf(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MAG_FILTER, GLES20.GL_LINEAR);
        GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_WRAP_S, GLES20.GL_CLAMP_TO_EDGE);
        GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_WRAP_T, GLES20.GL_CLAMP_TO_EDGE);
        //GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MIN_FILTER, GLES20.GL_NEAREST);
        //GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MAG_FILTER, GLES20.GL_NEAREST);

        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, mTextures[1]);

        GLES20.glTexParameterf(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MIN_FILTER, GLES20.GL_LINEAR);
        GLES20.glTexParameterf(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MAG_FILTER, GLES20.GL_LINEAR);
        GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_WRAP_S, GLES20.GL_CLAMP_TO_EDGE);
        GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_WRAP_T, GLES20.GL_CLAMP_TO_EDGE);
        //GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MIN_FILTER, GLES20.GL_NEAREST);
        //GLES20.glTexParameteri(GLES20.GL_TEXTURE_2D, GLES20.GL_TEXTURE_MAG_FILTER, GLES20.GL_NEAREST);

        camOrigTarg = initSmallPtGPU(1, 128, "kernels/rendering_kernel.cl", texW, texH, "scenes/obj-model.txt", Environment.getExternalStorageDirectory() + "/" + context.getPackageName(), context.getAssets());

        Log.d("VRApp3Renderer", "End of onSurfaceCreated");
    }

    public void onDrawFrame(HeadTransform headTransform, Eye leftEye, Eye rightEye) {
        Log.i("VRApp3Renderer", "Start of onDrawFrame");

        float camOrg[] = new float[3];
        float camTar[] = new float[3];

        camOrg[0] = camOrigTarg[0];camOrg[1] = camOrigTarg[1];camOrg[2] = camOrigTarg[2];
        camTar[0] = camOrigTarg[3];camTar[1] = camOrigTarg[4];camTar[2] = camOrigTarg[5];

        camOrg[0] -= 1;
        reinitCamera(camOrg[0], camOrg[1], camOrg[2], camTar[0], camTar[1], camTar[2]);

        int[] arr_pixels = updateRendering();
        if (arr_pixels == null)
        {
            Log.e("VRApp3Renderer", "null is returned!!!");
            return;
        }

        Bitmap leftBmp = Bitmap.createBitmap(arr_pixels, texW, texH, Bitmap.Config.ARGB_8888);

        camOrg[0] += 2;
        reinitCamera(camOrg[0], camOrg[1], camOrg[2], camTar[0], camTar[1], camTar[2]);

        arr_pixels = updateRendering();
        if (arr_pixels == null)
        {
            Log.e("VRApp3Renderer", "null is returned!!!");
            return;
        }

        Bitmap rightBmp = Bitmap.createBitmap(arr_pixels, texW, texH, Bitmap.Config.ARGB_8888);

        if (bFile) {
            saveToFile(leftBmp, rightBmp);
        }

        GLES20.glClear(GLES20.GL_DEPTH_BUFFER_BIT | GLES20.GL_COLOR_BUFFER_BIT);
        GLES20.glUseProgram(mProgramHandle);

        mScreenPosition.position( 0 );
        GLES20.glVertexAttribPointer(mPositionHandle, mPositionDataSize, GLES20.GL_FLOAT, false, 0, mScreenPosition);
        GLES20.glEnableVertexAttribArray(mPositionHandle);

        // Pass in the color information
        mTextureCoordinate.position( 0 );
        GLES20.glVertexAttribPointer(mTexCoordinateHandle, mTexCoordinateDataSize, GLES20.GL_FLOAT, false, 0, mTextureCoordinate);
        GLES20.glEnableVertexAttribArray(mTexCoordinateHandle);

        GLES20.glActiveTexture(GLES20.GL_TEXTURE0);

        // Draw a left bitmap
        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, mTextures[0]);
        GLES20.glUniform1i(mTextureUniformHandle, 0);

        GLUtils.texImage2D(GLES20.GL_TEXTURE_2D, 0, leftBmp, 0);
        GLES20.glDrawArrays(GLES20.GL_TRIANGLES, 0, 6);

        // Draw a right bitmap
        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, mTextures[1]);
        GLES20.glUniform1i(mTextureUniformHandle, 0);

        GLUtils.texImage2D(GLES20.GL_TEXTURE_2D, 0, rightBmp, 0);
        GLES20.glDrawArrays(GLES20.GL_TRIANGLES, 0, 6);

        // Finalize a drawing
        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, 0);
        GLES20.glUseProgram(0);

        Log.i("VRApp3Renderer", "End of onDrawFrame, "+String.valueOf(inv));
        inv++;
    }

    public void onDrawEye(Eye eye) {
        Log.i("VRApp3Renderer", "Start of onDrawEye");

        float camOrg[] = new float[3];
        float camTar[] = new float[3];

        camOrg[0] = camOrigTarg[0];camOrg[1] = camOrigTarg[1];camOrg[2] = camOrigTarg[2];
        camTar[0] = camOrigTarg[3];camTar[1] = camOrigTarg[4];camTar[2] = camOrigTarg[5];

        if (eye.getType() == Eye.Type.LEFT) {
            camOrg[0] -= 1;
        }
        else if (eye.getType() == Eye.Type.RIGHT) {
            camOrg[0] += 1;
        }

        reinitCamera(camOrg[0], camOrg[1], camOrg[2], camTar[0], camTar[1], camTar[2]);
        int[] arr_pixels = updateRendering();

        if (arr_pixels == null)
        {
            Log.e("VRApp3Renderer", "null is returned!!!");
            return;
        }

        Bitmap bitmap = Bitmap.createBitmap(arr_pixels, texW, texH, Bitmap.Config.ARGB_8888);

        if (bFile) {
            if (eye.getType() == Eye.Type.LEFT) saveToFile(bitmap, null);
            if (eye.getType() == Eye.Type.RIGHT) saveToFile(null, bitmap);
        }

        GLES20.glClear(GLES20.GL_DEPTH_BUFFER_BIT | GLES20.GL_COLOR_BUFFER_BIT);
        GLES20.glUseProgram(mProgramHandle);

        mScreenPosition.position( 0 );
        GLES20.glVertexAttribPointer(mPositionHandle, mPositionDataSize, GLES20.GL_FLOAT, false, 0, mScreenPosition);
        GLES20.glEnableVertexAttribArray(mPositionHandle);

        // Pass in the color information
        mTextureCoordinate.position( 0 );
        GLES20.glVertexAttribPointer(mTexCoordinateHandle, mTexCoordinateDataSize, GLES20.GL_FLOAT, false, 0, mTextureCoordinate);
        GLES20.glEnableVertexAttribArray(mTexCoordinateHandle);

        GLES20.glActiveTexture(GLES20.GL_TEXTURE0);
        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, mTextures[0]);
        GLES20.glUniform1i(mTextureUniformHandle, 0);

        GLUtils.texImage2D(GLES20.GL_TEXTURE_2D, 0, bitmap, 0);
        GLES20.glDrawArrays(GLES20.GL_TRIANGLES, 0, 6);

        // Finalize a drawing
        GLES20.glBindTexture(GLES20.GL_TEXTURE_2D, 0);
        GLES20.glUseProgram(0);

        Log.i("VRApp3Renderer", "End of onDrawEye, "+String.valueOf(inv));
        inv++;
    }

    @Override
    public void onSurfaceChanged(int width, int height) {
        Log.e("VRApp3Renderer", "Start of onSurfaceChanged");

        try {
            GLES20.glViewport(0, 0, width, height);//specifies transformation from normalized device coordinates to window coordinates
        } catch (Exception e) {
            // TODO Auto-generated catch block
            Log.d("Exception at the onSurfaceChanged",e.getMessage());
        }

        Log.e("VRApp3Renderer", "End of onSurfaceChanged");
    }

    @Override
    public void onFinishFrame(Viewport viewport) {
        finishRendering();
    }

    @Override
    public void onRendererShutdown() {

    }

    public void onNewFrame(HeadTransform headTransform) {
        float orig[] = new float[3];
        float targ[] = new float[3];

        headTransform.getHeadView(orig, 0);
        headTransform.getForwardVector(targ, 0);

        reinitCamera(orig[0], orig[1], orig[2], targ[0], targ[1], targ[2]);
    }
}
