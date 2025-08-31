#include "UI.h"


FPSCounter::FPSCounter() {
	frames = 0;
	currentFrame = 0;
	fps = 0.0;
	lastTime = std::chrono::steady_clock::now();
}

void FPSCounter::step() {
	frames++;
	auto now = std::chrono::steady_clock::now();
	double elapsed = std::chrono::duration<double>(now - lastTime).count();
	if (elapsed >= 1.0) {
		fps = frames / elapsed;
		frames = 0;
		lastTime = now;
	}
	currentFrame++;
}

double FPSCounter::getFPS() const { return fps; }
double FPSCounter::getFrameTime() const { return 1.0 / fps; }
int FPSCounter::getCurrentFrame() const { return currentFrame; }




Window::Window(unsigned int height, unsigned int width) {
    SCR_HEIGHT = height;
    SCR_WIDTH = width;

    double maxHeightWidth = std::max(height, width);
    windowScaleFactor = Vec3(height / maxHeightWidth, width / maxHeightWidth, 1);  // Scale factor that sacle x (-1.2->1.2) and y (-1->1) for example
    userScaleFactor = Vec3(1, 1, 1);
    sceneOrigin = Vec3(0, 0, 0);
    window = nullptr;
    uniform_varColour = NULL;

}


Window::~Window() {
    terminate();
}








bool Window::terminate() {

    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
    std::cout << "[INFO] VBO/VAO Buffers Successfully Destroyed" << std::endl;

    glfwTerminate();
    std::cout << "[INFO] GLFW Successfully Terminated" << std::endl;

    return 1;
}






bool Window::initGlfw() {

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif



    // glfw window creation
    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return 0;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    std::cout << "[INFO] GLFW Successfully Initiated" << std::endl;


    return 1;
}





bool Window::initGlad() {
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return 0;
    }
    std::cout << "[INFO] Glad Successfully Initiated" << std::endl;

    return 1;

}







bool Window::initFreetype() {
    FT_Library ft;
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
        return 0;
    }


    FT_Face face;
    if (FT_New_Face(ft, "./fonts/arial/ARIAL.TTF", 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
        return 0;
    }
    else {
        FT_Set_Pixel_Sizes(face, 0, 48);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        // load first 128 characters of ASCII set
        for (unsigned char c = 0; c < 128; c++)
        {
            if (FT_Load_Char(face, c, FT_LOAD_RENDER))
            {
                std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
                continue;
            }

            unsigned int texture;
            glGenTextures(1, &texture);
            glBindTexture(GL_TEXTURE_2D, texture);

            glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
            );

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            Character character = {
                texture,
                Vec3(face->glyph->bitmap.width, face->glyph->bitmap.rows),
                Vec3(face->glyph->bitmap_left, face->glyph->bitmap_top),
                static_cast<unsigned int>(face->glyph->advance.x)
            };
            Characters.insert(std::pair<char, Character>(c, character));
        }
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    // destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    std::cout << "[INFO] FreeType Successfully Initiated" << std::endl;

    return 1;

}







void Window::DrawText(std::string text, float x, float y, float scale, Vec3 col)
{
    fontShader.use();
    fontShader.setVec3("textColor", col);


    glActiveTexture(GL_TEXTURE0);

    std::cout << "draw text " << text << std::endl;

    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++)
    {
        Character ch = Characters[*c];

        float xpos = x + ch.Bearing.x * scale;
        float ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        float w = ch.Size.x * scale;
        float h = ch.Size.y * scale;


        float vertices[6][4] = {
            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos,     ypos,       0.0f, 1.0f },
            { xpos + w, ypos,       1.0f, 1.0f },

            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos + w, ypos + h,   1.0f, 0.0f }
        };



        /*glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);


        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glDrawArrays(GL_TRIANGLES, 0, 3);*/

        glBindTexture(GL_TEXTURE_2D, ch.TextureID);
        
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        glDrawArrays(GL_TRIANGLES, 0, 6);
        x += (ch.Advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    //glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}










void Window::SetScaleFactor(Vec3 scaleFactor) {
    userScaleFactor = scaleFactor.abs();
}

Vec3 Window::GetScaleFactor() {
    return userScaleFactor;;
}















void Window::DrawCircle(Circle circle, Vec3 col) {
    for (int i = 0; i < circle.triangles.size(); i++)
    {
        DrawTriangle(circle.triangles[i], col);
    }
}


void Window::DrawQuad(Quad quad, Vec3 col) {
    DrawTriangle(quad.triangles[0], col);
    DrawTriangle(quad.triangles[1], col);
}





void Window::DrawWireframeQuad(Quad quad, Vec3 col, float edgeThickness) {
    // Bottom edge (between x1 and x2)
    DrawQuad(
        Quad(
            quad.vertex1, quad.vertex2, Vec3(quad.vertex2.x, quad.vertex2.y + edgeThickness, quad.vertex2.z), Vec3(quad.vertex1.x, quad.vertex1.y + edgeThickness, quad.vertex1.z)
        ), col
    );
    // Top edge (between x4 and quad.vertex3)
    DrawQuad(
        Quad(
            Vec3(quad.vertex4.x, quad.vertex4.y - edgeThickness, quad.vertex4.z), Vec3(quad.vertex3.x, quad.vertex3.y - edgeThickness, quad.vertex3.z), quad.vertex3, quad.vertex4
        ), col
    );
    // Left edge (between quad.vertex1 and quad.vertex4)
    DrawQuad(
        Quad(
            quad.vertex1, Vec3(quad.vertex1.x + edgeThickness, quad.vertex1.y, quad.vertex1.z), Vec3(quad.vertex4.x + edgeThickness, quad.vertex4.y, quad.vertex4.z), quad.vertex4
        ), col
    );
    // Right edge (between quad.vertex2 and quad.vertex3)
    DrawQuad(
        Quad(
            Vec3(quad.vertex2.x - edgeThickness, quad.vertex2.y, quad.vertex2.z), quad.vertex2, quad.vertex3, Vec3(quad.vertex3.x - edgeThickness, quad.vertex3.y, quad.vertex3.z)
        ),
        col
    );
}





// Function to draw a triangle
void Window::DrawTriangle(Triangle tri, Vec3 col) {
    Vec3 combinedScaleFactor = Vec3(userScaleFactor.x * windowScaleFactor.x, userScaleFactor.y * windowScaleFactor.y, userScaleFactor.z * windowScaleFactor.z);

    float vertices[] = {
        tri.vertex1.x * combinedScaleFactor.x - sceneOrigin.x, tri.vertex1.y * combinedScaleFactor.y - sceneOrigin.y, tri.vertex1.z * combinedScaleFactor.z - sceneOrigin.z,
        tri.vertex2.x * combinedScaleFactor.x - sceneOrigin.x, tri.vertex2.y * combinedScaleFactor.y - sceneOrigin.y, tri.vertex2.z * combinedScaleFactor.z - sceneOrigin.z,
        tri.vertex3.x * combinedScaleFactor.x - sceneOrigin.x, tri.vertex3.y * combinedScaleFactor.y - sceneOrigin.y, tri.vertex3.z * combinedScaleFactor.z - sceneOrigin.z,
    };

    solidShader.use();
    solidShader.setVec3("colour", col);
    

    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glDrawArrays(GL_TRIANGLES, 0, 3);
}



bool Window::OnUserCreate() { return 1; }
bool Window::OnUserUpdate() { return 1; }



















bool Window::Init()
{

    // Init glfw Window
    if (!initGlfw()) { return 0; }

    // Init Glad
    if (!initGlad()) { return 0; }

    // Init freetype
    if (!initFreetype()) { return 0; }


    // build and compile our shader program
    solidShader.compile();
    fontShader.compile();


    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);



    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, this);
    glfwSetScrollCallback(window, scrollCallback);


    
    if (wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}




    if (!OnUserCreate()) {
        return 0;
    }



    while (!glfwWindowShouldClose(window))
    {
        
        processInput(window);


        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        //DrawText("SAMPLE TEXT", 1.0f, 0.0f, 1.0f, Vec3(0.5, 0.8f, 0.2f));


        if (!OnUserUpdate()) {
            return 0;
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    glBindVertexArray(0);


    if (!terminate()) { return 0; };

    return 1;

}




















// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void Window::processInput(GLFWwindow * window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    glfwGetCursorPos(window, &mousePos.x, &mousePos.y);

}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void Window::framebuffer_size_callback(GLFWwindow * window, int width, int height)
{
    glViewport(0, 0, width, height);
}



Vec3 Window::GetScreenOrigin() {
    return sceneOrigin;
}
void Window::SetScreenOrigin(Vec3 origin) {
    sceneOrigin = origin;
}




Vec3 Window::GetScreenSize() {
    return Vec3(SCR_WIDTH, SCR_HEIGHT, 0);
}


Vec3 Window::GetMousePosition() {
    return mousePos;
}



void Window::SetTitle(std::string title) {
	glfwSetWindowTitle(window, title.c_str());
}


void Window::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    Window* _this = (Window*)glfwGetWindowUserPointer(window);
    _this->scrollDirection = yoffset;
}


bool Window::IsScrollUp() {
    if (scrollDirection > 0) {
        scrollDirection = 0;
        return 1;
    }
    return 0;
}


bool Window::IsScrollDown() {
    if (scrollDirection < 0) {
        scrollDirection = 0;
        return 1;
    }
    return 0;
}
