#ifndef INCLUDED_GRAPHICS
#define INCLUDED_GRAPHICS 1
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <pthread.h>
#include <unistd.h>
#include <iostream>
namespace graphics
{
  class Frame;
  void playwindow(Frame* f);
  class Frame : public sf::RenderWindow
  {
  protected:
    int wid,hei;
    bool autoclear;
    double framerate=60;
    std::vector<sf::View> views;
  public:
    Frame(int w,int h,const std::string& title="SFML window",bool clr=true) : sf::RenderWindow(sf::VideoMode(w, h), title)
    {
      wid=w; hei=h; autoclear=clr;
      views.push_back(sf::View(sf::FloatRect(0.f,0.f, (float)w,(float)h)));
      //this->setView(views[0]);
    }

    virtual void onClose() {this->close();}
    virtual void onKeyPress(const sf::Event& e) {std::cout << "Key pressed: "<<e.key.code<<"\n";}
    virtual void onKeyRelease(const sf::Event& e) {std::cout << "Key released: "<<e.key.code<<"\n";}
    virtual void onMousePress(const sf::Event& e) {std::cout << "Mouse pressed\n";}
    virtual void onMouseRelease(const sf::Event& e) {std::cout << "Mouse released\n";}
    virtual void onMouseMoved(const sf::Event& e) {}
    virtual void onDraw() {}
    inline sf::Thread& launchWindow()
    {
      sf::Thread* thread=new sf::Thread(&playwindow, this);
      thread->launch();
      return *thread;
    }

    int addView(float cx,float cy,float w,float h) //Returns the index of the view (Rectangle construction - top-left corner, w, h)
    {
      views.push_back(sf::View(sf::FloatRect(cx,cy,w,h)));
      return views.size()-1;
    }
    inline void switchView(int ind) {this->setView(views[ind]);} //Switches to the view at i'th index
    inline void restoreView() {this->setView(views[0]);} //Switches back to full-view
    inline const std::vector<sf::View>& getAllViews() const {return views;}


    inline void setFrameRate(double d) {framerate=d;}
    friend void playwindow(Frame*);
  };
  void playwindow(Frame* f)
  {
    while (f->isOpen())
    {
        sf::Event event;
        while (f->pollEvent(event))
        {
            if (event.type == sf::Event::Closed) f->onClose();
            else if (event.type == sf::Event::KeyPressed) f->onKeyPress(event);
            else if (event.type == sf::Event::KeyReleased) f->onKeyRelease(event);
            else if (event.type == sf::Event::MouseButtonPressed) f->onMousePress(event);
            else if (event.type == sf::Event::MouseButtonReleased) f->onMouseRelease(event);
            else if (event.type == sf::Event::MouseMoved) f->onMouseMoved(event);
        }
        if(f->autoclear) f->clear();
        f->onDraw();
        f->display();
        double sleeptime=1000000/f->framerate;
        usleep(sleeptime);
    }
  }

  namespace objects
  {
    namespace shapes
    {
      inline static sf::CircleShape getCircle(float rad,const sf::Color& fill=sf::Color::White,const sf::Color& border=sf::Color::Black,float outline=0,double pX=0,double pY=0)
      {
        sf::CircleShape shape(rad);
        shape.setRadius(rad);
        shape.setFillColor(fill);
        shape.setOutlineThickness(outline);
        shape.setOutlineColor(border);
        shape.setPosition(pX,pY);
        return shape;
      }
      inline static sf::RectangleShape getCircle(float w, float h,const sf::Color& fill=sf::Color::White,const sf::Color& border=sf::Color::Black,float outline=0,double pX=0,double pY=0)
      {
        sf::RectangleShape shape(sf::Vector2f(w,h));
        shape.setFillColor(fill);
        shape.setOutlineThickness(outline);
        shape.setOutlineColor(border);
        shape.setPosition(pX,pY);
        return shape;
      }
    }
  }

}
#endif
