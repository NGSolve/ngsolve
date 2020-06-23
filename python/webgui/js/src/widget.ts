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
    } , 0);
    this.model.on('change:value', this.data_changed, this);
  }
  data_changed() {
    let render_data = this.model.get("value");
    this.scene.updateRenderData(render_data);
  }
}

