from dash import html, dcc, Input, Output, callback
import dash_mantine_components as dmc
from dash_iconify import DashIconify

def create_toc_sidebar():
    """
    Creates a fixed sidebar with a table of contents.
    """

    toc = html.Div(
        [
            dmc.Group(
                mb="md",
                children=[
                    DashIconify(icon="tabler:list-search", width=18),
                    dmc.Text("Table of contents"),
                ]
            ),
            html.Div(id="toc-container"),
        ],
        style={
            "position": "sticky",
            "top": "1rem"
        }
    )
    return toc

def register_toc_callbacks():
    @callback(
        Output("toc-container", "children"),
        Input("page-content", "children")
    )
    def update_toc(content):
        headings = []
        
        def extract_headings(component, depth=0):
            if component is None:
                return
                
            if isinstance(component, (str, int, float)):
                return
                
            if isinstance(component, dict):
                component_type = component.get('type', '')
                props = component.get('props', {})
                
                if (isinstance(component_type, str) and 
                    component_type.lower() in ['h1', 'h2', 'h3'] or 
                    getattr(component_type, '__name__', '').lower() in ['h1', 'h2', 'h3']):
                    
                    heading_text = props.get('children', '')
                    heading_id = props.get('id', '')
                    
                    if isinstance(heading_text, (list, dict)):
                        heading_text = str(heading_text)
                    
                    if heading_text and heading_id:
                        level = int(str(component_type)[-1]) if isinstance(component_type, str) else int(component_type.__name__[-1])
                        headings.append({
                            'level': level,
                            'text': heading_text,
                            'id': heading_id
                        })
                
                children = props.get('children', [])
                if children is not None:
                    if isinstance(children, (list, tuple)):
                        for child in children:
                            extract_headings(child, depth + 1)
                    else:
                        extract_headings(children, depth + 1)
            
            elif isinstance(component, (list, tuple)):
                for item in component:
                    extract_headings(item, depth + 1)

        extract_headings(content)
        
    #  <Box<'a'>
    #   component="a"
    #   href={item.link}
    #   onClick={(event) => event.preventDefault()}
    #   key={item.label}
    #   className={cx(classes.link, { [classes.linkActive]: active === item.link })}
    #   style={{ paddingLeft: `calc(${item.order} * var(--mantine-spacing-md))` }}
    # >
    #   {item.label}
    # </Box>

        toc_items = []
        for heading in headings:
            indent = (heading['level']) * 20
            toc_items.append(
                html.A(
                    heading['text'],
                    href=f"#{heading['id']}",
                    className="toc-item",
                    style={"paddingLeft": f"{indent}px"}
                )
            )
        
        return toc_items