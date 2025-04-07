import { MantineProvider } from '@mantine/core';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Layout from './layouts/Layout';
import Home from './pages/Home';
import Monomers from './pages/Monomers';
import MonomerDetail from './pages/MonomerDetail';
import About from './pages/About';
import ApiDocs from './pages/ApiDocs';

function App() {
  return (
    <MantineProvider withGlobalStyles withNormalizeCSS>
      <Router>
        <Layout>
          <Routes>
            <Route path="/" element={<Home />} />
            <Route path="/monomers" element={<Monomers />} />
            <Route path="/monomers/:id" element={<MonomerDetail />} />
            <Route path="/about" element={<About />} />
            <Route path="/api" element={<ApiDocs />} />
          </Routes>
        </Layout>
      </Router>
    </MantineProvider>
  );
}

export default App;
