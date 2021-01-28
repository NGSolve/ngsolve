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
  container: any;

  render() {
    // show preview image, a text message on hover
    // load real render data on click and start webgui
    let files = this.model.get("value");
    const image = files['preview'];
    this.container = $(`
      <div class="webgui_container" style="width:100%">
          <img src="${image}" class="image">
          <div class="webgui_overlay webgui_tooltip">
              <span class="webgui_tooltiptext"> Click to load interactive WebGUI </span>
          </div>
      </div>`);
    let div = document.createElement( 'div' );
    this.container.click( el=> this.onClickImage(el) )
    this.container.appendTo(div);
    this.el.appendChild(div);
    this.model.on('change:value', this.data_changed, this);
  }
  onClickImage(el) {
      document.body.style.cursor = "wait";
          let files = this.model.get("value");
          $.get(files['render_data'], (render_data) => {
              this.container.remove();
              this.container = null;
              document.body.style.cursor = "";
              let pel = this.el.children[0];
              pel.innerHTML = "";
              let scene = new Scene();
              scene.init(pel, render_data);
          });
  }

  data_changed() {
    let render_data = this.model.get("value");
    this.scene.updateRenderData(render_data);
  }
}


