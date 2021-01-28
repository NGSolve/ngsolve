// Copyright (c) Matthias Hochsteger
// Distributed under the terms of the Modified BSD License.

import {
  DOMWidgetView
} from '@jupyter-widgets/base';

import {
  Scene
} from './scene';

export class NGSolveView extends DOMWidgetView {
  scene: Scene;

  render() {
    console.log("Render NGSView");
    let render_data = this.model.get("value");
    console.log("render data", render_data);
    this.scene = new Scene();
    let container = document.createElement( 'div' );
    container.setAttribute("style", "height: 50vw; width: 100vw;");
    this.el.appendChild(container);
    setTimeout(()=> {
      this.scene.init(container, render_data);
      this.scene.render();
    } , 0);
    this.model.on('change:value', this.data_changed, this);
  }
  data_changed() {
    let render_data = this.model.get("value");
    this.scene.updateRenderData(render_data);
  }
}

export class NGSolveDocuView extends DOMWidgetView {
  scene: Scene;

  render() {
    console.log("Render NGSDocuView");
    let files = this.model.get("value");
    console.log("data", files);
    let container = document.createElement( 'div' );
    container.setAttribute("style", "height: 50vw;");
    let img = document.createElement('img');
    img.setAttribute("src", files['preview']);
    img.setAttribute("style", "width: 100%");
    img.onclick = (el) => this.onClickImage(el);
    container.appendChild(img);
    this.el.appendChild(container);
    this.model.on('change:value', this.data_changed, this);
  }
  onClickImage(el) {
      console.log("clicked image", el);
      document.body.style.cursor = "wait";
          let files = this.model.get("value");
          $.get(files['render_data'], (render_data) => {
              document.body.style.cursor = "";
              let pel = this.el.children[0];
              pel.innerHTML = "";
              console.log("got render data", render_data);
              let scene = new Scene();
              scene.init(pel, render_data);
          });
  }

  data_changed() {
    let render_data = this.model.get("value");
    console.log("got new render data", render_data);
    this.scene.updateRenderData(render_data);
  }
}


