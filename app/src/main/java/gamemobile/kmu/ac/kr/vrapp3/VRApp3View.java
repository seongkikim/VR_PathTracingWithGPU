package gamemobile.kmu.ac.kr.vrapp3;

import android.content.Context;
import android.view.MotionEvent;
import android.view.View;

import com.google.vr.sdk.base.GvrView;
import com.google.vr.sdk.base.GvrViewerParams;

import static android.view.MotionEvent.ACTION_DOWN;
import static android.view.MotionEvent.ACTION_UP;

/**
 * Created by user on 2018-03-22.
 */
public class VRApp3View extends GvrView{
    int touchX, touchY;
    private VRApp3Renderer mRenderer;
    private float distLens;

    public VRApp3View(Context context) {
        super(context);

        // Create an OpenGL ES 1.0 context.
        setEGLContextClientVersion(2);

        distLens = getGvrViewerParams().getInterLensDistance();
        mRenderer = new VRApp3Renderer(context, distLens);

        setRenderer(mRenderer);
        setTransitionViewEnabled(false);

        // Enable Cardboard-trigger feedback with Daydream headsets. This is a simple way of supporting
        // Daydream controller input for basic interactions using the existing Cardboard trigger API.
        enableCardboardTriggerEmulation();

        setOnTouchListener(new OnTouchListener() {
            @Override
            public boolean onTouch(View v, MotionEvent event) {
                if (event.getAction() == ACTION_DOWN) {
                    touchX = (int) event.getX();
                    touchY = (int) event.getY();
                }
                else if (event.getAction() == ACTION_UP) {
                    //mRenderer.onTouch(touchX - (int) event.getX(), touchY - (int) event.getY());
                }

                return true;
            };
        });
    }
}
